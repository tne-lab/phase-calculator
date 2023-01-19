/*
------------------------------------------------------------------

This file is part of a plugin for the Open Ephys GUI
Copyright (C) 2017 Translational NeuroEngineering Laboratory, MGH

------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef ASIC_PHASE_CALCULATOR_H_INCLUDED
#define ASIC_PHASE_CALCULATOR_H_INCLUDED

/*

The Phase Calculator filters selected input channels to the given passband,
generates an estimate of the phase of each band-limited signal in degrees
using the Hilbert transform, and outputs this as a continuous stream.
There are also options to output the magnitude or imaginary component of the analytic signal.

The visualizer panel serves a complementary role. When given a continuous channel
(that is enabled for processing, i.e. selected in the param window) and an event
channel to monitor, it calculates a delayed but precise estimate of the phase of the
filtered continuous channel at the time of each event onset and adds it to the rose plot.
If the events represent stimulation times, for example, this allows you to monitor the
accuracy of phase-locked stimulation in real time.

@see GenericProcessor, ASICPhaseCalculatorEditor

*/

#include <ProcessorHeaders.h>
#include "HTransformers.h"    // Hilbert transformers & frequency bands
#include <OpenEphysFFTW.h> 
#include <mutex>
#include <queue>
#include <DspLib.h> 

namespace ASICPhaseCalculator {
    
    static const float passbandEps = 0.01F;      // default passband width if pushing lowCut down or highCut up to fix invalid range, and also the minimum for lowCut. (meant to be larger than actual minimum floating-point eps)
    static const int glitchLimit = 200;          // "glitch limit >> How many samples should we check for unwrapping & smoothing the Caclulated Phase  
    static const int visHilbertLengthMs = 1024;  // process length for real-time visualization ground truth hilbert transform,based on evaluation of phase error compared to offline processing
    static const int visMinDelayMs = 675;
    static const int visMaxDelayMs = 1000;

    struct circularArray;     // Implementation of Circular Array for Fast Filtering
    struct ChannelInfo;       // Holds Information about all Channels          
	struct ActiveChannelInfo; // Holds Information/Data about Active Channels only. 
	class DataStreamSettings; // Holds the Settings for each DATA STREAMS
	class Node;               // The main Processor Class 
	class VisBufferManager;
	class CalVisualizerThread;
	class CompPhaseBufferManager;
    
    
	enum OutputMode { PH = 1, MAG, PH_AND_MAG, IM, NUM_OUTPUTMODES };  // Output mode - corresponds to itemIDs on the ComboBox
	enum EisText { NTEXT = 1, YTEXT = 2 };                             //added by sumedh
	enum EisLAA { NLAA = 1, YLAA = 2 };
	


    //%%%%%%%%%%%%%%%%%%%%%%%% CIRCULAR ARRAY STRUCTURE DECLARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Circular Array is used for Convolution [Filtering], instead of shifting for each sample we just put the new sample and index in circle.
    struct circularArray {
        Array <double> data;                         // Juce Array for storing elements of Circular Array
        int len = data.size();                       // Stores the Length of Array
        int head = 0;                                // Stores the Head [Latest Sample] of Circular Array       
        double* dataPtr = data.getRawDataPointer();  // Pointer to the Array [indexable]

        circularArray();                             // Constructor
        void resize(int size);                       // Resizing the Array
        void update();                               // Updating member functions [~ resetting]
        double at(int i);                            // Acquiring the Previous ith Sample
        void put(double item);                       // Adding a new value value into the Circular Array
    private:
        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(circularArray);
    };


	//%%%%%%%%%%%%%%%%%%%% BUFFER MANAGER FOR PHASE COMPARISIONS/VISUALIZATION  %%%%%%%%%%%%%%%%%%%%%%
	// This Class is Responsible for Visualizations....
	class VisBufferManager {
	public:
		int buffLen;                              // the Length of the Buffer [1.7s]
		int startIdx;                             // int(Downsampled-SampingRate*0.5); 
		int nextIdx;                              // The index of the Last Sample in the buffer
		int trigIdx;                              // After 1.5 Seconds Data is loaded The Comparision Should be Performed
		Array <double> visHilbertASIC;            // Hilbert Transform from ASIC Plugin.
		Array <double> visBPFdata;                // BAND PASS FILTERED DATA used for Phase Calculation using FFT Library
		bool visThreadHasCopiedData;              // is set to TURE if the Phase Comp Plugin has Copied the Contents
		std::queue<juce::int64> visTsBuffer;      // holds stimulation timestamps until the delayed phase is ready to be calculated
		juce::int64 lastBuffEndSampTS;            // Timestamp for Last Sample of Recent Incoming buffer
		
		std::mutex mutx;

		/** Constructor */
		VisBufferManager(int downSample2Rate);

		/** Destructor */
		~VisBufferManager();

		/** Resets/Resizes the BufferSizes/Indexes and Clears Buffers
			(Used when Plugin Settings Updated -> Sampling Rate)   */
		void reset(int downSample2Rate);

		/** Resets the indexes and Clears the buffer*/
		void clearBuffer();

		/** Used to insert new data into this class Buffers (Member Variables)*/
		void instertData(double* hilbertASIC, double* bpfData, int length , juce::int64 endTS);

		/** Used to Copy this Class Variables into the targetObj pass as argument*/
		void copyBuffer(VisBufferManager* targetObj);

		/** Returns true if the Buffer is triggered or (nextIdx > trigIdx) */
		bool isBufferLengthTriggered() { return  nextIdx > trigIdx; }

	private:
		JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(VisBufferManager);
	};



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION COMP THREAD CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This Class is Responsible for COMPARING ASIC PHASE VALUES WITH PHASE VALUES ACQUIRED FROM FFT TRANSFROM:
	class CalVisualizerThread : public Thread {
	public:

		/** Class Constructor */
		CalVisualizerThread(Node* nodePtr, int downSample2Rate);

		/** Class Destructor */
		~CalVisualizerThread();

		/** Main Function of The Thead*/
		void run() override;

		/** Used to Start this Phase Calc Thread */
		void startThisThread();

		/** Used to Stop this Thread */
		void stopThisThread();

		/** The Compared Phase Values are stored in a Queue
			Which is Shared with Visualizer thread for Rose Plot */
		void clearCompPhaseBuffer();

		/** Try to Read from the output Queue [Compared Phase Values] of this thread
			It check if Queue is not busy [No thread is using it], to read the Queue */
		bool tryToReadVisPhasesfromCompBuffer(std::queue<double>& targetQueue);

	private:

		Node* node;                                   // A Pointer to main ASIC Class Obj
		bool StopThreadFlag;                          // If TURE, breaks the main loop of Thread [ends it]
		VisBufferManager* visThreadBuffer;           // Buffer in Which Contents are Copied TO
		FFTWTransformableArray* visHilbertBufferFFT;  // The FFT buffer, used for aquiring FFT Phase Signal
		CompPhaseBufferManager* compPhaseBuffer;      // Pointer to Class that manages the Compared Phase Queue Vals

		JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(CalVisualizerThread);
	};


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUSALIZATION BUFFER MANAGER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This Class is Responsible for Managing the Final Visualization Buffer that Consists of Phase Values
	class CompPhaseBufferManager {
	public:
		CompPhaseBufferManager()
			: visPhaseBuffer()
			, visPhaseBufferCS() {
		}

		void insertNewPhaseData(float* newPhaseCompValues, int length) {
			ScopedLock phaseBufferLock(visPhaseBufferCS);
			for (int i = 0; i < length; i++) {
				visPhaseBuffer.push(double(newPhaseCompValues[i]));
			}
		}

		void clearBuffer() {
			ScopedLock phaseBufferLock(visPhaseBufferCS);
			while (!visPhaseBuffer.empty()) {
				visPhaseBuffer.pop();
			}
		}

		/* Locks and Moves the elements of visPhaseBuffer Queue into the Target Queue*/
		bool tryToReadVisPhases(std::queue<double>& targetQueue) {

			ScopedTryLock lock(visPhaseBufferCS);
			if (!lock.isLocked()) {
				return false;
			}

			while (!visPhaseBuffer.empty()) {
				targetQueue.push(visPhaseBuffer.front());
				visPhaseBuffer.pop();
			}
			return true;
		}

	private:
		std::queue<double> visPhaseBuffer;
		CriticalSection visPhaseBufferCS;
	};




    //%%%%%%%%%%%%%%%%%%%%%%%% ACTIVE CHANNEL INFO CLASS DECLARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Active Channel Info Structure is used to Keep Record of the Active Channels 
    struct ActiveChannelInfo {

        ActiveChannelInfo(const ChannelInfo* cInfo);   // Constructor
        void update();                                 // Updates Active Channel Info Class Object Memeber
        void reset();                                  // reset to perform after end of Acquisition or UPDATE

        circularArray cirLpfStateMain;     // Holds Samples for IIR Filter [Step I >> Node::Process]
        circularArray cirLpfState;         // Holds Samples for FIR LPF [Step III >> Node::Process]
        circularArray cirBpfState;         // Holds Samples for FIR BPF [Step IV >> Node::Process]
        circularArray cirhtState;          // Holds Samples for HILBERT FILTER (size = 31)
        circularArray cirHtRealState;      // Holds Samples for Hilbert Trasnform Real Part.

        int newlpfstartingpoint;     // During Sampling we have to keep track of the Next Buffer Starting Point [Default: 1000]           
        int newprelpfstartingpoint;  // During Sampling we have to keep track of the Next Buffer Starting Point [Default: 3000]
        double lastCompPhs2;         // Keeping Track of The Prev Phase Value [Because It might be in the Last Buffer]
        double lastCompMag2;         // Keeping Track of The Prev Mangnitude Value
        double nextCompPhase;        // Keeping Track of Next ComputePhase [The Prev Phase Needs to be Conveted into the Next Phase at each Buffer]
        double nextCompMag;          // Keeping Track of The Next Magnitude to Compute 
        int   strideLoc;             // Keeping Track of The Stride Location. 
        double phaseStep;            // Keeping Track of The Phase Steps.
        double magStep;              // Keeping Track of Magnitude Steps.
        float lastPhase;             // last phase output, for glitch correction [Unwrapping and Smoothing]
        const ChannelInfo* chanInfo; // Will Hold the Address to Channel Information Structure

        // VISUALLIZATION BUFFERS

        int indexVisBuffer;                       // Index for Visualization Buffer [Each Time we have to append Samples]
        int visNsamples;                          // Length of Visualization buffers [# Samples for 1.5 seconds data]
        float lastPhaseValueVis;                  // Used for Unwrapping
        Array <double> visBPFdata;                // 1.5 Seconds Long LPF Buffer 
        Array <float> visHilbertBufferASIC;       // 1.5 Seconds Long ASIC Hilbert Visualization Buffer


    private:
        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ActiveChannelInfo);
    };


    //%%%%%%%%%%%%%%%%%%%%%%%% CHANNEL INFO STRUCT DECLARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Channel Info Class is used to Keep Record of All Channels 
    struct ChannelInfo {

        int chan;                 // index of channel among owner's data channels
        bool isActivated;         // Check if the Channel is Activated
        int dsFactor;             // 0 if sample rate is not a multiple of Hilbert::FS (in this case it cannot be activated.)
        float sampleRate;         // Channel Sampling Rate
        const Node* node;         // Data Stream

        void update();            // Updates 
        bool activate();          // returns false on failure. does nothing if already active.
        void deactivate();        // Deactivates the Channel
        bool isActive() const;    // Checks the Status of the Channel

        std::unique_ptr<ActiveChannelInfo> acInfo;    // info for ongoing phase calculation - null if non-active.
        ChannelInfo(Node* asicNode,  int i);          // Constructor. Requires Node Object and Channle Index Number

    private:
        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ChannelInfo);
    };

   
    
	class Node : public GenericProcessor
	{
	public:
		/** The class constructor, used to initialize any members. */
		Node();

		/** The class destructor, used to deallocate memory */
		~Node();

		/** If the processor has a custom editor, this method must be defined to instantiate it. */
		AudioProcessorEditor* createEditor() override;

		/** Called every time the settings of an upstream plugin are changed [DEFAULT] */
		void updateSettings() override;

		/***/
		void createEventChannels() override;

		/** Defines the functionality of the processor.
			The process method is called every time a new data buffer is available.
			Visualizer plugins typically use this method to send data to the canvas for display purposes */
		void process(AudioBuffer<float>& buffer) override;

		/** Handles events received by the processor
			Called automatically for each received event whenever checkForEvents() is called from
			the plugin's process() method */
			//void handleTTLEvent(TTLEventPtr event) override;

			/** Handles spikes received by the processor
				Called automatically for each received spike whenever checkForEvents(true) is called from
				the plugin's process() method */
		//void handleSpike(SpikePtr spike) override;

		/** Handles broadcast messages sent during acquisition
			Called automatically whenever a broadcast message is sent through the signal chain */
		//void handleBroadcastMessage(String message) override;

		/** Saving custom settings to XML. This method is not needed to save the state of
			Parameter objects */
			//void saveCustomParametersToXml(XmlElement* parentElement) override;

			/** Load custom settings from XML. This method is not needed to load the state of
				Parameter objects*/
				//void loadCustomParametersFromXml(XmlElement* parentElement) override;

				/** Retruns the Selected Band by the ASIC Processor,
					This Funtion is also called by the Editor to Populate the Band Value*/
		Band getBand() const;

		/** Returns the Ouput Mode Value*/
		OutputMode getOutputMode() const;

		/** Gets the DownSampling 1 Rate*/
		int getDownSample1Rate() const;

		/** Gets the DownSampling 2 Rate*/
		int getDownSample2Rate()  const;

		/** Gets the Value for isLAA */
		EisLAA getIsLAA() const;

		/** Gets the Value for isText*/
		EisText getIsText() const;

		/** Returns Index of the Visualization Channel*/
		int getVisContinuousChannelIndex() {
			return visContinuousChannel;
		}


		/** Gets the Status of EDITOR LOCK
			If returns true the the EDITOR SHOULD BE LOCKED*/
		bool getLockEditorInfo() const;
		void setLockEditor(bool newState);

		/** Sets (Updates) the Frequency Band
			The Band Number is Passed as Argument, {Band>> type:enum}	*/
		void setBand(Band newBand);

		/** Sets (Updates) the OUTPUT MODE of ASIC Plugin [enum]
			Modes: PHASE, MAGNITUDE, PHASE AND MAGNITUDE*/
		void setOutputMode(OutputMode newOutputMode);


		/** Includes the MAGNITUDE Channels If Ouputmode is Set to PH_AND_MAG */
		void updateMagnitudeChannels(bool updateDownwardStreams);



		/** Sets (Updates) FIRST STAGE SAMPLING RATE,
			DOWN SAMPLING 1 UPDATE METHOD */
		void setDownSample1Rate(int newSamplingRate);

		/** Sets (Updates) Second STAGE SAMPLING RATE,
			DOWN SAMPLING 2 UPDATE METHOD */
		void setDownSample2Rate(int newSamplingRate);

		/** Sets (Updates) The Selection for Liear Approximation
			Updates isLAA:  Use linear Approximation or not*/
		void setIsLAA(EisLAA newIsLAA);

		/** Sets (Updates) The Selection for TxT outputs
			Updates isText: Print intermidate outputs to Txt File or Not*/
		void setIsText(EisText newIsText);

		/** This method is Called Whenever the Channel Parameter is Changed
			Called by the GenericProcessor::setParameter())*/
		//void parameterValueChanged(Parameter* param) override;

		/** This Function is executed At the Start of the Acquisition.
			The Return Type should be Ture for this Function for Starting the Aquisition
			Here this Function is used to ENABLE the LOCK [Variable] to lock the Editor and Prevent PARAM UPDATE */
		//bool startAcquisition() override;


		/** This Function is executed At the END of the Acquisition
			The Return type should be ture for sucessfully stoping Acquisition*/
		//bool stopAcquisition() override;

		// Coefficients of LOW PASS FILTER!
		Array<double> transFun = Array<double>({ -0.00173272639840494,-0.00198748297232062,-0.00242153251168736,-0.00267113412739108,-0.00211925399627538,1.23339961810775e-18,0.00443805687582359,0.0117422605199169,0.0220815503798583,0.0351289228347335,0.0500324747539425,0.0654894059431507,0.0799158415661605,0.0916844398106024,0.0993861724102343,0.102066009823313,0.0993861724102343,0.0916844398106024,0.0799158415661605,0.0654894059431507,0.0500324747539425,0.0351289228347335,0.0220815503798583,0.0117422605199169,0.00443805687582359,1.23339961810775e-18,-0.00211925399627538,-0.00267113412739108,-0.00242153251168736,-0.00198748297232062,-0.00173272639840494 });
		static double FilterMATLAB(double input, circularArray& state, int nCoefs, const double* transF);

		// Use LAA to return phase angle in degrees
		double LAA(std::complex<double>);

		/** Circular distance between angles x and ref (in radians). The "cutoff" is the maximum
			possible positive output; greater differences will be wrapped to negative angles.
			For interpolation, and also RosePlot visualization.*/
		static double circDist(double x, double ref, double cutoff = 2 * Dsp::doublePi);

		// Do glitch unwrapping
		void unwrapBuffer(float* wp, int nSamples, float lastPhase);

		// Do start-of-buffer smoothing
		void smoothBuffer(float* wp, int nSamples, float lastPhase);


		/** This Function is called for Updating Continous Channels as More Channels are added into the Stream
			Note DataStream and Processor Keep Seperate Records for Constinous Channels and we update both*/
		void updateChannelsFor_PH_AND_MAG();


		// Retruns the ID for the Selected Stream 
		uint16 getSelectedStreamID();

		

		/** Returns array of ACTIVE channels that only includes inputs(not extra outputs)
			The Indexes of ACTIVE CHANNELS are Returned*/
		Array<int> getActiveInputs() const;

		/** Convenience method to call "update" on all active channel info structs in the map.
			Resets All of the Active Channels to Default/Initial Values */
		void updateAllActiveChannels();

		/** Enables an input channel to be processed. Returns false if this is prohibited due to
			the channel's sample rate. If the channel is already active, does nothing.
			Precondition: the channel passed in is a valid input channel index. */
		bool activateInputChannel(int chan);

		/** This Functions Deactivates the Input Channel
			The Channel index is Pass to this Functions for Deactivation*/
		void deactivateInputChannel(int chan);

		/** Clears Current Channel Array for the Specified [DataStream (DS)] stored in [DataStreamSettings]
			Then initializes [Channel Info] Array storing information for each Channel in the DS
			Whever Input DS Changes This function should be Called for setting up CH INFO for each CH in a DS*/
		void resetAllChInfos4DataStream();

		/** Renews All Active Channels [Resets + Updates Active Channels]
			When Band Is Changed This needs to be RUN to update Buffer Sizes
			Doesn't Recreate Active Channels, Just Renew/Rest Current Active Channels*/
		void renewAllActiveChannels();

		/** Recreate All Active Channels [Creates Active Info Classes for Active Channel]
			The Record of Active Channels is Stored as An Array of Ints in ASIC Class */
		void recreateAllActiveChannels(Array <int>, int totalChannels);


		int getInputChCount() const;
		void setInputChCount(int newCount);


		bool isGeneratesTimestamps() const override;
		int getNumSubProcessors() const override;
		float getSampleRate(int subProcessorIdx = 0) const override;
		float getBitVolts(int subProcessorIdx = 0) const override;

		bool enable() override;
		bool disable() override;


		// maps full source subprocessor IDs of incoming streams to indices of
		// corresponding subprocessors created here.
		HashMap<int, uint16> subProcessorMap;

		// Update subProcessorMap
		void updateSubProcessorMap();


		//-----------------MADE PUBLIC ---------------------------
		OwnedArray<ChannelInfo> channelInfo;
		Array <int> selActiveChannels;
		uint16 selectedStreamID;
		//--------------------------------------------------------

	
		//%%%%%%%%%%%%%%%%%%%%%%%%%% ASIC VISUALIZATION METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		// miscellaneous helper [Accessed by Canvas.Cpp File]
		int getFullSourceId(int chan);

		/** Changes/Sets the Visualization Continous Channel by setting [visCotinousChannel]
			TO VERIFY ->It also remove the Timestamps of the Previous Channel stored in [visTsBuffer] Queue */
		void setVisContChan(int newChan);

		/** Changes/Sets the Visualization Channel by setting [visEventChannel]
			TO VERIFY->It also remove the Timestamps of the Previous Channel stored in [visTsBuffer] Queue */
		void setVisEventChan(int newChan);


		/** Check the Visualization Timestamp Queue, Clearing the ones that have expired
			(Late to Calc Phase), While Caculate phase of any that are ready. sdbEndTs = timestamp 1
			Past the end of the current buffer.*/
		void calcVisPhases(ActiveChannelInfo* acInfo, juce::int64 sdbEndTs);


		/** Reads from the visPhaseBuffer if it can acquire a TryLock. Returns True if Successful*/
		bool tryToReadVisPhases(std::queue<double>& other);

		/** Respond to Incoming events if a stimEventChannel is selected.*/
		//void handleTTLEvent(TTLEventPtr event) override;


		/** Get Visualization ASIC Buffer Pointer */
		VisBufferManager* getVisBufferPtrASIC() { return visASICBuffer; }

		/** Get Pointer to Visuliazation Phase Compare Thread*/
		CalVisualizerThread* getCalVisualizerThreadPtr() {
			return calVisualizerThreadPtr.get();
		}

		/** Start the Visualization Thread*/
		void startCalVisualizerThread() {
			calVisualizerThreadPtr->startThisThread();
		}
		
		/** Stop the Visualization Thread*/
		void stopCalVisualizerThread() {
			calVisualizerThreadPtr->stopThisThread();
		}

		/** For Passing an Even to Crossing Detector Optimization*/
		EventChannel* getVisPhaseChannel() { 
			return visPhaseChannel; 
		}

		/** Allows the other thread to access Node Class Private Function*/
		void addEventThread(juce::int64 LastSampleTS, double eventData) {
			juce::int64 eventTs  = LastSampleTS - getNumSamples(visContinuousChannel);
			BinaryEventPtr event = BinaryEvent::createBinaryEvent(visPhaseChannel, eventTs, &eventData, sizeof(double));
			addEvent(visPhaseChannel, event, 0);
		}

	private:


		// Allow responding to stim events if a stimEventChannel is selected.
		void handleEvent(const EventChannel* eventInfo, const MidiMessage& event, int samplePosition = 0) override;


		int inputChannelCount;   // Keeps a Record for the input Channel Count
		Band band;               // frequency band (determines which Hilbert transformer to use) [enum]
		OutputMode outputMode;   // Output Mode (determines which Mode we are using i.e PHASE, MAGNITUDE, PHASE_AND_MAGNITUDE, etc)[enum]
		int dsIscale;            // Down Sampling  I Rate 
		int dsIIscale;           // Down Sampling II Rate
		EisLAA isLAA;            // isLAA: Is LAA Approximation Enabled
		EisText isText;          // isText: Is Printing to TXT files Enabled
		bool lockEditor;         // Locks the Editor while the Acquisition is Reunning.

		// --------------------------------------------------------------------------
		//uint16 selectedStreamID; // The ID of the Selected Stream [Made Public for OLD OE VERSION]
		//Array <int> selActiveChannels; [Made Public for OLD OE VERSION]
		//---------------------------------------------------------------------------

		/** For MAG_AND_PH, Additional Channels are introduced, But Processor and Datastream Keep a Record of the Seperately.
			But when the Dropbox of Output Mode is Changed, The gereral listener (ParameterValueChaned) is executed. So The
			MAG Channels are added by both the CHANNELS and the UpdateSettings Function which result in extra copies of datastream channels
			WE NEED TO UPDATED IT ONLY ONES*/
		bool magChannelsNeedsUpdate;


		// for Canvas
		int visEventChannel;      // Event Channel to Watch for Phase Plotting (-1 = No Channel)
		int visContinuousChannel; // Channel to calculate Phases from at receieved stim event times.
		std::unique_ptr <CalVisualizerThread> calVisualizerThreadPtr;
		VisBufferManager* visASICBuffer;

		// event channel to send visualized phases over
		EventChannel* visPhaseChannel;

	};


}






#endif // PHASE_CALCULATOR_H_INCLUDED
