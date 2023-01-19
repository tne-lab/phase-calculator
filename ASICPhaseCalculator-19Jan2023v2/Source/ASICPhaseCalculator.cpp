/*
------------------------------------------------------------------

This file is part of a plugin for the Open Ephys GUI
Copyright (C) 2022 Translational NeuroEngineering Laboratory, MGH

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

#include <chrono>    // For Calculating Timings [Speed]
#include <cfloat>    // DBL_MAX
#include <cmath>     // sqrt
#include <cstring>   // memcpy, memmove
#include <windows.h>
#include <string.h>
#include <fstream>
#include "ASICPhaseCalculator.h"
#include "ASICPhaseCalculatorEditor.h"

namespace ASICPhaseCalculator {



    //%%%%%%%%%%%%%%%%%%%%%%%% CIRCULAR ARRAY CLASS FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Circular Array is used for Convolution [Filtering], instead of shifting for each sample we just put the new sample and index in circle.
    
    circularArray::circularArray() {                              // Constructor 
        update();                                                 // Calls the Update Function for initialization of the member varibles.
    }        
    
    void circularArray::update() {                                // To Update Class Variables, This is done because when the CircularArray [CircularArray::data] is resized it needs an Update 
        len = data.size();  head = 0;                             // Update the Length and set the Head at the Start.
        FloatVectorOperations::clear(data.begin(), data.size());  // Clears the Data Vector
        dataPtr = data.getRawDataPointer();                       // Pointer to the Vector Data
    }

    void circularArray::resize(int size) {                        // To Resize the Length of Circular Array
        data.resize(size);                                        // Resizing the Size of Array
        update();                                                 // Calling the Update Function as other member variables requires an update.
    }   
    
    double circularArray::at(int i) {                             // Getting the ith  Value
        if (head - i < 0)                                         // Going Backwards, if we Pass 0th indx Start from Rigth Most End
            return dataPtr[head - i + len];                       // as head - i < 0 so getting values from Right Most End towards Left
        else
            return dataPtr[head - i];                             // Return Values if there are still some of the left
    }

    void circularArray::put(double item) {
        if (head == len - 1)                                      // If Last Time We Had Updated The END, Then Start from the Beginning
            head = 0;                                             // Set Head to 0
        else
            ++head;                                               // Increment the Head
        dataPtr[head] = item;                                     // Add Recent Value into the Head
    }

    //%%%%%%%%%%%%%%%%%%%%%%% BUFFER MANAGER FOR VISUALIZATION MEMEBER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VisBufferManager::VisBufferManager(int downSample2Rate) {
        reset(downSample2Rate);
    }

    VisBufferManager::~VisBufferManager() {
        clearBuffer();//clears the queue
        visHilbertASIC.clear(); 
        visBPFdata.clear();
    }

    void VisBufferManager::reset(int downSample2Rate) {
        buffLen = int(1.7 * downSample2Rate);
        startIdx = int(0.5 * downSample2Rate);
        trigIdx = int(1.5 * downSample2Rate);
        nextIdx = startIdx;
        visThreadHasCopiedData = false;
        visHilbertASIC.resize(buffLen);
        visBPFdata.resize(buffLen);
        clearBuffer();
    }

    void VisBufferManager::clearBuffer() {
        double* visBPFdataPtr = visBPFdata.getRawDataPointer();
        double* visHilbertASICPtr = visHilbertASIC.getRawDataPointer();
        for (int i = 0; i < buffLen; i++) {
            visBPFdataPtr[i] = 0;
            visHilbertASICPtr[i] = 0;
        }
        while (!visTsBuffer.empty()) visTsBuffer.pop();// Also Empty the TimeStamp Queue  
    }

    void VisBufferManager::copyBuffer(VisBufferManager* targetObj) {
        std::lock_guard lock(mutx);
        targetObj->buffLen = buffLen;
        targetObj->startIdx = startIdx;
        targetObj->trigIdx = trigIdx;
        targetObj->nextIdx = nextIdx;

        targetObj->lastBuffEndSampTS = lastBuffEndSampTS;

        double* sourceVisBPFdataPtr = visBPFdata.getRawDataPointer();
        double* sourceVisHilbertASICPtr = visHilbertASIC.getRawDataPointer();

        double* destVisBPFdataPtr = targetObj->visBPFdata.getRawDataPointer();
        double* destVisHilbertASICPtr = targetObj->visHilbertASIC.getRawDataPointer();

        // Copying The Buffers
        for (int i = 0; i < nextIdx; i++) {
            destVisBPFdataPtr[i] = sourceVisBPFdataPtr[i];
            destVisHilbertASICPtr[i] = sourceVisHilbertASICPtr[i];
        }

        // Also Copying the TTL Event Time Stamps
        while (!visTsBuffer.empty()) {
            targetObj->visTsBuffer.push(visTsBuffer.front());
            visTsBuffer.pop();
        }

        visThreadHasCopiedData = true;
    }


    void VisBufferManager::instertData(double* hilbertASIC, double* bpfData, int length , juce::int64 endTS) {
        std::lock_guard lock(mutx);

        double* visBPFdataPtr = visBPFdata.getRawDataPointer();
        double* visHilbertASICPtr = visHilbertASIC.getRawDataPointer();

        lastBuffEndSampTS = endTS;

        // IF Buffer is Filled [Upto TrigIdx] and The Visualization Phase Comparision Thread
        // Has Copied The Buffer. Then we Shift The Buffer by 500 Samples to the Left.
        if (nextIdx > trigIdx && visThreadHasCopiedData) {
            for (int i = 0; i < startIdx; i++) {
                visBPFdataPtr[i] = visBPFdataPtr[startIdx + i];
                visHilbertASICPtr[i] = visHilbertASICPtr[startIdx + i];
            }
            for (int i = 2 * startIdx; i < nextIdx; i++) {
                visBPFdataPtr[i - startIdx] = visBPFdataPtr[i];
                visHilbertASICPtr[i - startIdx] = visHilbertASICPtr[i];
            }
            nextIdx = nextIdx - startIdx;
            for (int i = 0; i < length; i++) {
                visBPFdataPtr[nextIdx + i] = bpfData[i];
                visHilbertASICPtr[nextIdx + i] = hilbertASIC[i];
            }
            nextIdx = nextIdx + length;
            visThreadHasCopiedData = false;
        }

        // Other Wise we Fill the Buffer with the New Data from the ASIC Plugin.
        else {
            for (int i = 0; i < length; i++) {

                // If the Phase Comparision Thread Doesn't works Properly then the Buffer will Overflow will Occur
                if (nextIdx + i >= buffLen) {    // If the Buffer is filled and the other Thread hasnt copied yet
                    std::cout<<" [BUFFER OVERFLOW!] The Thread for Phase Calculation is Not Working Properly"<<std::endl;
                    jassertfalse;
                }

                visBPFdataPtr[nextIdx + i] = bpfData[i];
                visHilbertASICPtr[nextIdx + i] = hilbertASIC[i];
            }
            nextIdx = nextIdx + length;
        }
    }



    //%%%%%%%%%%%%%%%%%%%%%%%% ACTIVE CHANNEL INFO STRUCT MEMBER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Active Channel Info Class is used to Keep Record of the Active Channels 
    ActiveChannelInfo::ActiveChannelInfo(const ChannelInfo* cInfo)
        : chanInfo(cInfo) {                  // Constructor, The Object cInfo [ChannelInfo cInfo] is passed to Active info
        update();
    }

    void ActiveChannelInfo::update() {                                    // Used to Update the Class Object Variables
        const Node* p = chanInfo->node;
        cirLpfStateMain.resize(31);                                   // Resizing the IIR Low Pass Filter (31 Samples Filtering)
        cirLpfState.resize(lowpassfilt::delay[p->getBand()]);         // Resizing the FIR Low Pass Filter 
        cirBpfState.resize(bandpassfilt::delay[p->getBand()]);        // Resizing the FIR band Pass Filter
        cirhtState.resize(Hilbert::delay[p->getBand()]);              // Resiging the Hilbert Transform Filter
        cirHtRealState.resize(int(Hilbert::delay[p->getBand()] / 2)); // Delay to introduce in the Real Part [Orignal Input] of Hilbert Transform Signal 

        //indexVisBuffer = 0;                                // Index Location for Storing Visualization Buffer Data.
        //visNsamples = int(p->getDownSample2Rate() * 1.5);  // # of Samples for 1.5 Seconds Worth of Data
        //visBPFdata.resize(visNsamples);                    // Resizing for 1.5 Seconds Worth of Data [BPF downsampled data (Default:1000S/S)]
        //visHilbertBuffer.resize(visNsamples);              // Resizing for 1.5 Seconds Worth of Data [FFT HILBERT data for Comparision (Default:1000S/S)]
        //visHilbertBufferASIC.resize(visNsamples);          //

        reset();
    }

    void ActiveChannelInfo::reset() {
        cirLpfStateMain.update();   // ReNew/Reset the Circular Arrays
        cirLpfState.update();
        cirBpfState.update();
        cirhtState.update();
        cirHtRealState.update();

        newlpfstartingpoint = 0;     // Records Starting Position of Nextbuffer for Sampling [default 1000]
        newprelpfstartingpoint = 0;  // Records Starting Position of Nextbuffer for Sampling [default 3000]

        lastPhase = 0;      // Will Hold Last Sample Value of Previous Buffer
        lastCompPhs2 = 0;   // [Used in Interpolation] Will Hold the last Downsampled Phase
        lastCompMag2 = 0;   // [Used in Interpolation] Will Hold the last Downsampled Magnitude
        nextCompPhase = 0;  // [Used in Interpolation] Will Hold the next Downsampled Phase
        nextCompMag = 0;    // [Used in Interpolation] Will Hold the next Downsampled Magnitude
        strideLoc = 1;      // Starting start from [1 to Stride]. (Ignoring Index 0)
        phaseStep = 0;      // Starting Phase Step is 0
        magStep = 0;        // Starting Magitude Step is 0
    }

    //%%%%%%%%%%%%%%%%%%%%%%%% CHANNEL INFO STRUCT DECLARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Channel Info Structure is used to Keep Record of the ALL Channels 

    ChannelInfo::ChannelInfo(Node* asicNode, int i)     // Constructor
        : chan(i)                                       // Reqires a Channel Number
        , sampleRate(0)
        , dsFactor(0)
        , acInfo(nullptr)
        , node(asicNode)
        , isActivated(false)
    {
        update();
    }


    void ChannelInfo::update() {                                        // Update ChannelInfo Object Variables

        const DataChannel* contChannel = node->getDataChannel(chan);   // Acquire Pointer to the Data Channel Object
        if (contChannel == nullptr) {
            jassertfalse;
            return;
        }

        sampleRate = contChannel->getSampleRate();               // Acquire Sampling Rate of the Channel
        float fsMult = sampleRate / Hilbert::fs;
        float fsMultRound = std::round(fsMult);
        if (std::abs(fsMult - fsMultRound) < FLT_EPSILON) {
            dsFactor = int(fsMultRound);                         // can be active - sample rate is multiple of Hilbert Fs
            if (isActivated) {                                   // Updates/Resets Active Channels
                acInfo->update();
            }
        }
        else {
            dsFactor = 0;
            deactivate();                                        // this channel can no longer be active.
        }
    }


    bool ChannelInfo::activate() {                     // Activates the Channel if Not Active
        if (!isActivated && dsFactor != 0) {           // IF Channel Not Active
            acInfo.reset(new ActiveChannelInfo(this)); // Create an instance of ActiveChannelInfo and Store a Pointer to it.
            isActivated = true;
            std::cout << "[Channel Info] CH # " << chan << " ACTIVATED" << std::endl;
        }
        return isActivated;                            // Return weather the Channel is active or not
    }


    void ChannelInfo::deactivate() {                   // For deactivation Set the ActiveChannelInfo Pointer to NULL
        if (isActivated) {
            acInfo.reset();
            isActivated = false;
            std::cout << std::endl<<"[Channel Info] CH # " << chan << " DEACTIVATED" << std::endl;
        }
    }


    bool ChannelInfo::isActive() const {               // Check if the Channel is Active 
        return isActivated;
    }


    //%%%%%%%%%%%%%%%%%%%%%%%% DataStreamSettings CLASS FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // 
    // Retruns Array of Active Channels Indexes [Array of Int(Channel Indexes)]
    Array<int> Node::getActiveInputs() const {
        Array<int> activeInputs;                             // The indexes of Active Channels will be Stored in this Array
        for (auto chanInfo : channelInfo) {                  // Looping Through All Channels to find which ones are Active
            if (chanInfo->isActive()) {
                activeInputs.add(chanInfo->chan);            // If the Channel is Active Append its index to ativeInputs
            }
        }
        return activeInputs;
    }

    // Convenience method to call "update" on all active channel info structs in the map.
    void Node::updateAllActiveChannels() {
        for (int ai : getActiveInputs()) {
            jassert(channelInfo[ai] && channelInfo[ai]->isActive());
            channelInfo[ai]->acInfo->update();
        }
    }

    void Node::renewAllActiveChannels() {
        for (int ai : getActiveInputs()) {
            jassert(channelInfo[ai] && channelInfo[ai]->isActive());
            channelInfo[ai]->acInfo->reset();
        }
    }

    // Activates the Selected Input Channel, [Channel index is Passed to the Function]
    bool Node::activateInputChannel(int chan) {
        if (chan < 0 || chan >= channelInfo.size()) {
            jassertfalse;
            return false;
        }
        //jassert(!channelInfo[chan]->isActive());      
        //LOGD(" ");
        std::cout << "[Activate Function] CH # " << chan << " Active Before: " << channelInfo[chan]->isActive()<< std::endl;
        //LOGD("[Activate Function] CH #", chan, " Active Before: ", channelInfo[chan]->isActive());
        return channelInfo[chan]->activate();
    }


    // Deactivates the Selected Channel [Channel index is Passed for Deactivation]
    void Node::deactivateInputChannel(int chan) {
        if (chan < 0 || chan >= channelInfo.size()) {
            jassertfalse;
            return;
        }
        //jassert(channelInfo.getUnchecked(chan)->isActive());
        //LOGD(" ");
        std::cout << "[Deactivate Function] CH #" << chan << " Active Before: " << channelInfo[chan]->isActive() << std::endl;
        //LOGD("[Deactivate Function] CH #", chan, " Active Before: ", channelInfo[chan]->isActive());
        channelInfo.getUnchecked(chan)->deactivate();
    }

    // Clears and Reinstantiates Channel INFO array for a Specific Data Stream passed as Argument
    void Node::resetAllChInfos4DataStream() {
        channelInfo.clear();                                             // Clearing Channel Info Array
        for (int i = 0; i < inputChannelCount; ++i) {
            ChannelInfo* addChanInfo = new ChannelInfo(this, i);         // Adding New Channels in Channel Info Array
            channelInfo.add(addChanInfo);
        }
    }

    void Node::recreateAllActiveChannels(Array <int> selectedChannels, int totalChannels) {
        for (int i = 0; i < totalChannels; i++) {
            if (selectedChannels.contains(i)) {
                activateInputChannel(i);
            }
            else {
                deactivateInputChannel(i);
            }
        }
    }

    int Node::getInputChCount() const {
        return inputChannelCount;
    }

    void Node::setInputChCount(int newCount) {
        inputChannelCount = newCount;
    }



    //%%%%%%%%%%%%%%%%%%%%%%%% NODE CLASS DECLARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // NODE Is the Main Processor Class of the Plugin
  
    Node::Node() : GenericProcessor("ASIC") {

        setProcessorType(PROCESSOR_TYPE_FILTER);

        setOutputMode(PH);         // Setting the ASIC Processor OUTPUT for PHASE MODE
        setIsLAA(NLAA);            // Setting the ASIC Processor Not to USE LAA First
        setIsText(NTEXT);          // Setting the ASIC Processor Not to Print Intermidiate Outputs to TxT file
        setDownSample1Rate(3000);  // Setting the Down Sampling 1 Rate to 3000
        setDownSample2Rate(1000);  // Setting the Down Sampling 2 Rate to 1000
        lockEditor = false;        // Unlocking the Editor Initially 

        selActiveChannels.clear(); // Clears the Array [will store active Channels]
        selectedStreamID = 0;      // Store the Value for the Selected Stream ID 

        band = THETA;              // Setting BAND, [set Function not used] because BAND Change needs [Channle Info and Active Info]
                                   // Updates called in the function, However, initally they are NULL so BAND is initiazlied Directly.

        magChannelsNeedsUpdate = false; // Used to prevent Channels Update multiple Times 

        calVisualizerThreadPtr = std::make_unique < CalVisualizerThread >(this, getDownSample2Rate());
        visASICBuffer = new VisBufferManager(getDownSample2Rate());

        visContinuousChannel = -1;
        visEventChannel = -1; 
    }



    Node::~Node() {
        delete visASICBuffer;
        delete editor;
    }                                   

   
    AudioProcessorEditor* Node::createEditor(){       // Create a new instance of the Editor
        editor = new Editor(this);
        return editor;
    }

    
    
    void Node::createEventChannels(){                 // Used for creating Event Channel
        const DataChannel* visChannel = getDataChannel(visContinuousChannel);
        if (!visChannel){
            visPhaseChannel = nullptr;
            return;
        }

        float sampleRate   = visChannel->getSampleRate();
        EventChannel* chan = new EventChannel(EventChannel::DOUBLE_ARRAY, 1, 1, sampleRate, this);
        chan->setName(chan->getName() + ": PC visualized phase (deg.)");
        chan->setDescription("The accurate phase in degrees of each visualized event");
        chan->setIdentifier("phasecalc.visphase");

        // metadata storing source data channel
        MetaDataDescriptor sourceChanDesc(MetaDataDescriptor::UINT16, 3, "Source Channel",
            "Index at its source, Source processor ID and Sub Processor index of the channel that triggers this event",
            "source.channel.identifier.full");
        MetaDataValue sourceChanVal(sourceChanDesc);
        uint16 sourceInfo[3];
        sourceInfo[0] = visChannel->getSourceIndex();
        sourceInfo[1] = visChannel->getSourceNodeID();
        sourceInfo[2] = visChannel->getSubProcessorIdx();
        sourceChanVal.setValue(static_cast<const uint16*>(sourceInfo));
        chan->addMetaData(sourceChanDesc, sourceChanVal);

        visPhaseChannel = eventChannelArray.add(chan);
    }
    

   

    // %%%%%%%%%%%%%%% MAIN PROCESS FUNCTION CALLED EACH TIME THE BUFFER IS UPDATED %%%%%%%%%%%%%%%
    void Node::process(AudioSampleBuffer& buffer)
    {
        //std::ofstream myfileinput;
        //myfileinput.open("TIMINGS.txt", std::ofstream::out | std::ofstream::app);
        //auto start = std::chrono::high_resolution_clock::now();
        
        
        // handle subprocessors, if any
        HashMap<int, uint16>::Iterator subProcIt(subProcessorMap);
        while (subProcIt.next()){
            uint32 fullSourceID = uint32(subProcIt.getKey());
            int subProcessor = subProcIt.getValue();
            uint64 sourceTimestamp = getSourceTimestamp(fullSourceID);
            uint32 sourceSamples = getNumSourceSamples(fullSourceID);
            setTimestampAndSamples(sourceTimestamp, sourceSamples, subProcessor);
        }
        
        // check for events to visualize
        bool hasCanvas = static_cast<Editor*>(getEditor())->canvas != nullptr;
        if (hasCanvas && visEventChannel > -1){
            checkForEvents();
        }
       
        // iterate over active input channels
        Array<int> activeChans = getActiveInputs();     // activeChans >> [Indexes of Channels that are Active]
        int numActiveChans = activeChans.size();        // numActiveChans >> Number of Active Channels
		for (int ac = 0; ac < numActiveChans; ++ac){    

			ChannelInfo* chanInfo = channelInfo[activeChans[ac]];  // Allocating Memory for Array of channelInfo[Structure] for the Active Channels only
            ActiveChannelInfo* acInfo = chanInfo->acInfo.get();          // ChannelInfo[Structure] includes ActiveChannelInfo[Structure], so, acInfo is a pointer to Active Channel Structure Instance
            
			int chan = chanInfo->chan;              // chan >> Channel INDX
			int nSamples = getNumSamples(chan);     // nSamples >> # Samples
			if (nSamples == 0){                     // If there is no samples in the Buffer, Nothing TO DO
				continue;
			}

			int downsampleA = dsIscale;                           // dsIscale >> Node Class Variable [Holds First Downsampling Rate]
			int downsampleB = dsIIscale;                          // dsIIscale >> Node Class Variable [Holds Second Downsampling Rate]
            std::ofstream myfileinput;                            // Text File Handler Object [Used for printing the outputs into text files]
            float* const wpInput = buffer.getWritePointer(chan);  // Pointer to the Start of the Input Buffer [Holds the samples]

            // isText == 2, Print the Input Data onto a txt File
			if (isText == 2) {          
				myfileinput.open("inputdata.txt", std::ofstream::out | std::ofstream::app);
				for (int i = 0; i < nSamples; i++){
					myfileinput << wpInput[i] << "," << std::endl;
				}
				myfileinput.close();
			}

            // Step I: IIR Low pass Filter over entire data 
            double* newLPF = new  double[nSamples];     // Allocating memory for Storing Filtered Output Signal
            for (int i = 0; i < nSamples; i++) {
                newLPF[i] = FilterMATLAB(wpInput[i], acInfo->cirLpfStateMain, 31, transFun.begin());
                int j = i;
            }

            // Step II: Down Sampling I. Sampling the input data to [3000 samples/second : Default]
            std::vector<double> newdslpf;                                               // Vector for storing down Sampled A Signal [First Down Sampling]
            std::vector<int> idxRecDwnSamp_A;                                           // Recording Index of the Input Signal Selected for Down Sampling [Will be used for interpolation later] 
            int preDSample = int(ceil(chanInfo->sampleRate / dsIscale));                // Calculating the down sampling Stride
            int newprelpfstartingpoint = acInfo->newprelpfstartingpoint;                // Aquiring the Stating Point of Sampling [As data comes in Small buffers, 
            int prelpfjumps = 0;                                                        // so the Sampling Starting Point Depents of where we left Sampling in the last buffer]
            for (int i = newprelpfstartingpoint; i < nSamples; i = i + preDSample) {
                idxRecDwnSamp_A.push_back(i);                                           // Storing index of Sampled Value
                newdslpf.push_back(newLPF[i]);                                          // Storing Sampled Signal                        
                prelpfjumps = i;
            }
            if (prelpfjumps < nSamples) {                                               // If the final Sampled value [at index = prelpfjumps] is not the last value of the buffer
                acInfo->newprelpfstartingpoint = (preDSample + prelpfjumps) - nSamples; // then sample the next buffer after preDSample Samples as if the two buffers are cancatenated. 
            }

            // if isText == 2, Print the Down Sampled A data onto the TXT file
            if (isText == 2) {
                myfileinput.open("LPFdownsampleI.txt", std::ofstream::out | std::ofstream::app);
                for (int i = 0; i < newdslpf.size(); i++) {
                    myfileinput << newdslpf[i] << std::endl;
                }
                myfileinput.close();
            }

            // Step III: FIR low pass filter from MATLAB UPDATE
            double* newlpfData = new  double[newdslpf.size()];       // Allocating memory for Storing FIR Filtered Signal
            for (int lpfIndex = 0; lpfIndex < newdslpf.size(); ++lpfIndex) {
                newlpfData[lpfIndex] = FilterMATLAB(newdslpf[lpfIndex], acInfo->cirLpfState, lowpassfilt::delay[band], lowpassfilt::transformer[band].begin());
            }

            // Step IV: Down Sampling II. Sampling the Sampled Data Again [default 1000 samples/second]
            int new_sr_size = dsIscale / dsIIscale;                 // Caculating the Down Sampling B Stride [Sampling Jump Value]
            std::vector<double> newdsLpfData;                       // Vector for Storing Down Sampled B values.[Second Down Sampling]
            std::vector<int> idxRecDwnSamp_B;                       // Recording Index Locations of Down Sampling B values with respect to the input Signal.[Used for interpolation]
            int newlpfstartingpoint = acInfo->newlpfstartingpoint;  // Acquired the Starting Point of Buffer Sampling [Based on the previous buffer right most unsampled values]
            int lpfjumps = 0;
            for (int i = newlpfstartingpoint; i < newdslpf.size(); i = i + new_sr_size) {
                idxRecDwnSamp_B.push_back(idxRecDwnSamp_A[i]);                             // Appending DownSampled B values
                newdsLpfData.push_back(newlpfData[i]);                                     // Appending Indexes for DownSampled B values
                lpfjumps = i;
            }
            if (lpfjumps < newdslpf.size()) {                                              // Calculating Sampling Start index of the Next Buffer, [incase some samples are left at the end]
                acInfo->newlpfstartingpoint = (new_sr_size + lpfjumps) - newdslpf.size();
            }

            // if isText == 2, Print the Down Sampled B data onto the TXT file
            if (isText == 2) {
                myfileinput.open("LPFdownsampleII.txt", std::ofstream::out | std::ofstream::app);
                for (int i = 0; i < newdsLpfData.size(); i++) {
                    myfileinput << newdsLpfData[i] << std::endl;
                }
                myfileinput.close();
            }

            // Step V: FIR Band pass filtering on the downsampled data
            double* newdsBPFdata = new  double[newdsLpfData.size()];                  // Allocating memory for Storing FIR BP Filtered Signal
            for (int bpfIndex = 0; bpfIndex < newdsLpfData.size(); ++bpfIndex) {
                newdsBPFdata[bpfIndex] = FilterMATLAB(newdsLpfData[bpfIndex], acInfo->cirBpfState, bandpassfilt::delay[band], bandpassfilt::transformer[band].begin());
            }

            // Step VI: calculate Hilbert transform
            int htOutputSamps = newdsLpfData.size();
            Array<std::complex<double>> htOutput;
            Array <double> visHilbertASIC;  //<< Visualization BUFFER
            htOutput.resize(newdsLpfData.size());
            for (int hilbertIndex = 0; hilbertIndex < htOutputSamps; ++hilbertIndex) {
                double ic = FilterMATLAB(newdsBPFdata[hilbertIndex], acInfo->cirhtState, Hilbert::delay[band], Hilbert::transformer[band].begin());
                double rc = acInfo->cirHtRealState.at(6);                 // Get the Past 7th Sample = 6th indx 
                acInfo->cirHtRealState.put(newdsBPFdata[hilbertIndex]);   // Put the current Sample in Circular Array
                htOutput.set(hilbertIndex, std::complex<double>(rc, ic));
                visHilbertASIC.add(arg(htOutput[hilbertIndex]));
            }

            if (isText == 2) {
                myfileinput.open("hilbertTransform-Real-Imag-Phase.txt", std::ofstream::out | std::ofstream::app);
                for (int i = 0; i < htOutput.size(); i++)
                    myfileinput << std::real(htOutput[i]) << "," << std::imag(htOutput[i]) << "," << std::arg(htOutput[i]) << std::endl;
                myfileinput.close();
            }

            // Step VII: INTERPOLATION
            int stride2 = preDSample * new_sr_size;       // Stride for interpolation
            double* outputPhase = new  double[nSamples];  // To Record Up-Sampled Phase
            double* outputMagni = new  double[nSamples];  // To Record Up-Sampled Magnitude
            double prevCompPhase = acInfo->lastCompPhs2;  // Previous Computed Phase [From the Last Buffer]
            double prevCompMag = acInfo->lastCompMag2;    // Previous Computed Magnitude [From the Last Buffer]
            double nextCompPhs2 = acInfo->nextCompPhase;  // Latest Computed Phase of the Last Buffer
            double nextCompMag2 = acInfo->nextCompMag;    // Latest Computed Magnitude of the Last Buffer
            int    indexStride = acInfo->strideLoc;       // Index Location of the Stride we left in the Last Buffer
            double phaseStep2 = acInfo->phaseStep;        // Phase Step from Previous Buffer
            double magStep2 = acInfo->magStep;            // Magnitude Step from Previous Buffer
            bool  needPhs2 = outputMode != MAG;
            bool  needMag2 = outputMode != PH;


            for (int i = 0, j = 0; i < nSamples; i++, indexStride++) {
                if (i == idxRecDwnSamp_B[j]) {
                    indexStride = 1;
                    if (needPhs2) {
                        prevCompPhase = nextCompPhs2;                                                        // Update the Previous Phase 
                        nextCompPhs2 = (isLAA == NLAA) ? std::arg(htOutput[j]) : LAA(htOutput[j]);          // Get the Next PhaseSample 
                        phaseStep2 = circDist(nextCompPhs2, prevCompPhase, Dsp::doublePi) / stride2;         // Small Phase steps between 2 Samples
                    }
                    if (needMag2) {
                        prevCompMag = nextCompMag2;
                        nextCompMag2 = std::abs(htOutput[j]);
                        magStep2 = (nextCompMag2 - prevCompMag) / stride2;     // Default 30 Steps (or 30 + 1 Pieces if First and Last Includede)
                    }
                    if (j < idxRecDwnSamp_B.size() - 1)
                        ++j;  // Incrementing for the Next Sample for all Downsampled Input Signal
                }

                if (needPhs2)
                    outputPhase[i] = circDist(prevCompPhase + phaseStep2 * indexStride, 0, Dsp::doublePi);
                if (needMag2)
                    outputMagni[i] = prevCompMag + magStep2 * indexStride;

            }
            acInfo->lastCompPhs2 = prevCompPhase;   // Updating Values for the next incoming buffer
            acInfo->lastCompMag2 = prevCompMag;
            acInfo->nextCompPhase = nextCompPhs2;
            acInfo->nextCompMag = nextCompMag2;
            acInfo->strideLoc = indexStride;
            acInfo->phaseStep = phaseStep2;
            acInfo->magStep = magStep2;


            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // NOTE THE LOW PASS FILTER IS ADDING GAIN IN THE SIGNAL So, I am using HACK
            float addGain = 0.04;// ADDING CONST GAIN FOR MAGNITUDE AND IMAGINARY OUTPUT
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



            // Step VIII: MAPPING THE PHASE and MAGNITUDE VALUES ONTO THE OUTPUT BUFFER
            float* wpOut = buffer.getWritePointer(chan);
            float* wpOut2;
            if (outputMode == PH_AND_MAG) {                                                           // second output channel

                //int totChanlsGenericProcessor = getTotalContinuousChannels();
                //int NChannels = getDataStream(selectedStreamID)->getChannelCount();
                int x = buffer.getNumChannels();
                int y = buffer.getNumSamples();
                int z = getNumOutputChannels();
                int w = getNumInputChannels();
                int v = getNumInputs();
                int u = getNumOutputs();


                wpOut2 = buffer.getWritePointer(getInputChCount() + ac);  // NoN Mag Channels + Active Channels

            }


            for (int i = 0; i < nSamples; i++) {
                switch (outputMode) {
                case MAG:
                    wpOut[i] = float(outputMagni[i]) * addGain;
                    break;
                case PH_AND_MAG:
                    wpOut2[i] = float(outputMagni[i]) * addGain;
                    // fall through
                case PH:
                    wpOut[i] = outputPhase[i] * 57.2958; // Converted into Degress from Radians
                    break;
                case IM:
                    wpOut[i] = float(outputMagni[i] * std::sin(outputPhase[i])) * addGain;
                    break;
                }
            }


            // Step IX unwrapping / smoothing
            if (outputMode == PH || outputMode == PH_AND_MAG) {
                unwrapBuffer(wpOut, nSamples, acInfo->lastPhase);
                smoothBuffer(wpOut, nSamples, acInfo->lastPhase);
                acInfo->lastPhase = wpOut[nSamples - 1];
            }


            if (isText == 2) {
                myfileinput.open("Output.txt", std::ofstream::out | std::ofstream::app);
                for (int i = 0; i < nSamples; i++)
                    myfileinput << wpOut[i] << std::endl;
                myfileinput.close();
            }


            /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
            /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
            // TO BE REMOVED

            /*
            myfileinput.open("TimeStamps.txt", std::ofstream::out | std::ofstream::app);
            for (int i = 0; i < nSamples; i++) {
                myfileinput << wpInput[i] << "," << getNumSamplesInBlock(selectedStreamID) << "," << getFirstSampleNumberForBlock(selectedStreamID) << "," << getFirstTimestampForBlock(selectedStreamID)<<std::endl;
            }
            myfileinput.close();
            */


            /*

            double* hilBufferPtr = acInfo->visHilbertBuffer.getRealPointer();
            int startIndex = acInfo->indexVisBuffer;
            int i;
            for (i = 0; i < newdsLpfData.size() && i+ startIndex < acInfo->visNsamples; i++) {
                hilBufferPtr[startIndex+i] = newdsBPFdata[i];
                acInfo->visBPFdata.insert(startIndex+i, newdsBPFdata[i]);
                acInfo->visHilbertBufferASIC.insert(startIndex + i, arg(htOutput[i]));
                //
            }
            acInfo->indexVisBuffer = startIndex + i;


            if (startIndex == acInfo->visNsamples) {

                unwrapBuffer(acInfo->visHilbertBufferASIC.getRawDataPointer(), acInfo->visNsamples, 0);
                smoothBuffer(acInfo->visHilbertBufferASIC.getRawDataPointer(), acInfo->visNsamples, 0);



                acInfo->visHilbertBuffer.hilbert();
                myfileinput.open("VisualizationTrials.txt", std::ofstream::out | std::ofstream::app);
                for (int i = 0; i < acInfo->visNsamples; i++) {
                    myfileinput << i << "," << arg(acInfo->visHilbertBuffer.getAsComplex(i)) * 57.2958 << "," << acInfo->visBPFdata[i] << ","<< acInfo->visHilbertBufferASIC[i] * 57.2958 << std::endl;
                }
                myfileinput.close();
                acInfo->indexVisBuffer = 0;
            }


            */






            /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
            /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


            // If this is the same channel as that of visualization channel then add visualization data
            if (hasCanvas && chan == visContinuousChannel) {
                // And Perform Visualization Stuff
                juce::int64 enTS = getTimestamp(chan) + getNumSamples(chan); // Last Sample Time Stamp
                visASICBuffer->instertData(visHilbertASIC.getRawDataPointer(), newdsBPFdata, newdsLpfData.size(), enTS);
            }

            







            // Deleting the Space Allocated on HEAP
            delete[] newLPF;
            delete[] newlpfData;
            delete[] newdsBPFdata;
            delete[] outputPhase;
            delete[] outputMagni;
        }
    }

    double Node::FilterMATLAB(double input, circularArray& state, int nCoefs, const double* transf) {
        state.put(input);
        double retValue = 0.0;
        for (int kCoef = 0; kCoef < nCoefs; ++kCoef) {
            retValue = retValue + state.at(kCoef) * transf[kCoef];
        }
        return retValue;
    }

    
    double Node::circDist(double x, double ref, double cutoff) {
        static const double twoPi = 2 * Dsp::doublePi;
        double xMod = std::fmod(x - ref, twoPi);
        double xPos = (xMod < 0 ? xMod + twoPi : xMod);
        return (xPos > cutoff ? xPos - twoPi : xPos);
    }

    // %%%%%%%%%%%%% UWRAPPING THE PHASE JUMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    void Node::unwrapBuffer(float* wp, int nSamples, float lastPhase) {                     // This buffer unwraps the sudden Phase Jumps in the Output Data
        for (int startInd = 0; startInd < nSamples - 1; startInd++) {                        // Indexing through all of the samples in the buffer one by one [Possible Start of the Jump]
            float diff = wp[startInd] - (startInd == 0 ? lastPhase : wp[startInd - 1]);     // Finding the Phase Difference between the Current Sample and the Previous Sample
            if (abs(diff) > 180) {                                                           // Search forward for a jump in the opposite direction
                int endInd;
                int maxInd;
                if (diff < 0) {                                                              // for Downward Jumps, unwrap if there's a jump back up within glitchLimit samples
                    endInd = -1;                                                            // endInd >> index where the Jump ends [Set if to -1, we will find the end location of the jump later]
                    maxInd = jmin(startInd + glitchLimit, nSamples - 1);                    // maxInd >> The maximum number of samples to look after a Sudden change is recorded.  
                }
                else {                                                                       // for upward jumps, default to unwrapping until the end of the buffer,
                    endInd = nSamples;                                                      // but stop if there's a jump back down sooner.
                    maxInd = nSamples - 1;
                }
                for (int currInd = startInd + 1; currInd <= maxInd; currInd++) {             // Finding the end index of the JUMP
                    float diff2 = wp[currInd] - wp[currInd - 1];
                    if (abs(diff2) > 180 && ((diff > 0) != (diff2 > 0))) {                   // If the end Index is found save it and beak the loop [we don't want to look further]
                        endInd = currInd;
                        break;
                    }
                }
                for (int i = startInd; i < endInd; i++) {                                    // Unrapping from [startInd, endInd)
                    wp[i] -= 360 * (diff / abs(diff));
                }
                if (endInd > -1) {                                                          // skip to the end of this unwrapped section
                    startInd = endInd;
                }
            }
        }
    }

    // %%%%%%%%%%%%% SMOOTHING AT THE START OF BUFFER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    void Node::smoothBuffer(float* wp, int nSamples, float lastPhase) {
        int actualGL = jmin(glitchLimit, nSamples - 1);
        float diff = wp[0] - lastPhase;
        if (diff < 0 && diff > -180) {                                             // identify whether signal exceeds last sample of the previous buffer within glitchLimit samples.
            int endIndex = -1;
            for (int i = 1; i <= actualGL; i++) {
                if (wp[i] > lastPhase) {
                    endIndex = i;
                    break;
                }
                else if (wp[i] - wp[i - 1] < -180 && (wp[i] + 360) > lastPhase) { // corner case where signal wraps before it exceeds lastSample
                    wp[i] += 360;
                    endIndex = i;
                    break;
                }
            }
            if (endIndex != -1) {                                                 // interpolate points from buffer start to endIndex
                float slope = (wp[endIndex] - lastPhase) / (endIndex + 1);
                for (int i = 0; i < endIndex; i++) {
                    wp[i] = lastPhase + (i + 1) * slope;
                }
            }
        }
    }



    void Node::updateSettings() {
        int numInputs = getNumInputs();         // returns total Number of Inputs to PLUGIN  
        int prevNumInputs = channelInfo.size(); // Previous # of Inputs [Lenght of channelInfo]

        int nToRemove = jmax(prevNumInputs - numInputs, 0);  // If we have removed some inputs
        channelInfo.removeLast(nToRemove);
        setInputChCount(numInputs);


        resetAllChInfos4DataStream();    // Resets All Channels for This Data Stream
        recreateAllActiveChannels(selActiveChannels, getInputChCount());

        updateSubProcessorMap();
        // Update PH and MAG Channels
        //updateChannelsFor_PH_AND_MAG();

        visASICBuffer->reset(getDownSample2Rate());
        calVisualizerThreadPtr->clearCompPhaseBuffer();
    }


    Band Node::getBand() const { return band; }                   // Returns the Recent band By ASIC Processor
    OutputMode Node::getOutputMode() const { return outputMode; } // Returns the Recent Ouput Mode
    int Node::getDownSample1Rate() const { return dsIscale; }     // Returns the Recent Down Sampling 1 Rate
    int Node::getDownSample2Rate() const { return dsIIscale; }    // Returns the Recent Down Sampling 1 Rate
    EisLAA Node::getIsLAA()  const { return isLAA; }              // Returns the Recent isLAA value [Approximation Check]
    EisText Node::getIsText() const { return isText; }            // Returns the Recent isText value [Intermidate Output to Txt file Check]
    bool Node::getLockEditorInfo() const { return lockEditor; }   // Returns the EIDTOR LOCK STATUS
    void Node::setLockEditor(bool state) { lockEditor = state; }

    void Node::setBand(Band newBand) {
        if (newBand == band) { return; }           // Return if the newBand is the Same as Old Band
        if (newBand < 0 || newBand >= NUM_BANDS) { // Check if the Band [Emumerator Value] is valid
            jassertfalse; return;
        }
        band = newBand;                            // Update the Band
        std::cout << "[ASIC]->[FREQ BAND] SET TO " << band << std::endl;

        // >>>>>>>>>>>>>>>>  NOTE: STREAM COUNT SHOULD BE ONE FOR NOW <<<<<<<<<<<<<<<<<<<<<<<
        renewAllActiveChannels();
    }

    void Node::setOutputMode(OutputMode newOutputMode) {
        //if (outputMode == newOutputMode) { return; }
        if (newOutputMode < 1 || newOutputMode >= NUM_OUTPUTMODES) {    // Check if the Band [Emumerator Value] is valid
            jassertfalse; return;
        }

        // In PH_AND_MAG settings we add extra channels which needs to be updated.
        if (outputMode == PH_AND_MAG || newOutputMode == PH_AND_MAG) {
            outputMode = newOutputMode;
            CoreServices::updateSignalChain(editor.get());
            std::cout << "[ASIC]->[OUTPUT MODE] SET TO " << outputMode << std::endl;
        }
        else {
            outputMode = newOutputMode;
            std::cout << "[ASIC]->[OUTPUT MODE] SET TO " << outputMode << std::endl;
        }

    }


    void Node::setDownSample1Rate(int newSamplingRate) {
        dsIscale = newSamplingRate;
        std::cout << "[ASIC]->[DOWN SAMPLE 1] SET TO " << dsIscale << std::endl;
    }
    void Node::setDownSample2Rate(int newSamplingRate) {
        dsIIscale = newSamplingRate;
        std::cout << "[ASIC]->[DOWN SAMPLE 2] SET TO " << dsIIscale << std::endl;
    }
    void Node::setIsLAA(EisLAA newIsLAA) {
        isLAA = newIsLAA;
        std::cout << "[ASIC]->[is LAA] SET TO " << isLAA << std::endl;
    }
    void Node::setIsText(EisText newIsText) {
        isText = newIsText;
        std::cout << "[ASIC]->[is TXT] SET TO " << isText << std::endl;
    }



    
    bool Node::isGeneratesTimestamps() const
    {
        return true;
    }

    int Node::getNumSubProcessors() const
    {   
        return subProcessorMap.size();
    }

    float Node::getSampleRate(int subProcessorIdx) const
    {
        jassert(subProcessorIdx < getNumSubProcessors());
        int chan = getDataChannelIndex(0, getNodeId(), subProcessorIdx);
        return getDataChannel(chan)->getSampleRate();
    }

    float Node::getBitVolts(int subProcessorIdx) const
    {
        jassert(subProcessorIdx < getNumSubProcessors());
        int chan = getDataChannelIndex(0, getNodeId(), subProcessorIdx);
        return getDataChannel(chan)->getBitVolts();
    }


    double Node::LAA(std::complex<double> c) {
        double q;   // Hold phase before quantization

        // Determine phase based on octant. See LAA alg details.
        if (std::abs(c.real()) >= std::abs(c.imag())) {

            if (c.real() >= 0)
                q = (1. / 8.) * (c.imag() / c.real());       // octant 1 and 8
            else if (c.imag() >= 0)
                q = .5 + (1. / 8.) * (c.imag() / c.real());  // octant 4
            else
                q = -.5 + (1. / 8.) * (c.imag() / c.real()); // octant 5
        }
        else {
            if (c.imag() >= 0)
                q = 0.25 - (1. / 8.) * (c.real() / c.imag()); // octant 2 and 3
            else
                q = -.25 - (1. / 8.) * (c.real() / c.imag()); // octant 6 and 7
        }

        // Do quantization on phase (based on bit percision). We are testing 8 bits.
        double lsb = 1. / pow(2., 8.);
        double ans = floor(q / lsb) * lsb * 3.14 / 0.5;
        if (std::isnan(ans))
            ans = 0.0;
        return ans;
    }
 

    bool Node::enable() {                    // starts thread when acquisition begins
        if (isEnabled) {
            Editor* editor = static_cast<Editor*>(getEditor());  // have to manually enable editor, I guess...
            editor->enable();
        }
        return isEnabled;
    }


    bool Node::disable() {
        Editor* editor = static_cast<Editor*>(getEditor());
        editor->disable();

        for (auto chanInfo : channelInfo) {      // RESET states of active inputs
            if (chanInfo->isActive()) {
                chanInfo->acInfo->reset();
            }
        }

        // clear timestamp and phase queues
        //while (!visTsBuffer.empty())
        //{
        //    visTsBuffer.pop();
        //}

        while (!visASICBuffer->visTsBuffer.empty()) {
            visASICBuffer->visTsBuffer.pop();
        }

        calVisualizerThreadPtr->clearCompPhaseBuffer();

        return true;
    }

    // Handles the Event of Incoming Sham Pulses from Comparator Plugin
    void Node::handleEvent(const EventChannel* eventInfo, const MidiMessage& event, int samplePosition) {
        if (visEventChannel < 0) {
            return;
        }

        if (Event::getEventType(event) == EventChannel::TTL) {
            TTLEventPtr ttl = TTLEvent::deserializeFromMessage(event, eventInfo);
            if (ttl->getChannel() == visEventChannel && ttl->getState())  {
                
                juce::int64 ts = ttl->getTimestamp();                              // add timestamp to the queue for visualization
                
                //jassert(visTsBuffer.empty() || visTsBuffer.back() <= ts);
                //visTsBuffer.push(ts);

                jassert (visASICBuffer->visTsBuffer.empty() || visASICBuffer->visTsBuffer.back() <= ts);
                visASICBuffer->visTsBuffer.push(ts);
            }
        }
    }

    
    void Node::updateSubProcessorMap() {    // fill map according to selected channels, and remove outdated entries.

        if (outputMode != PH_AND_MAG) {
            subProcessorMap.clear();
            return;
        }

        uint16 maxUsedIdx = 0;
        SortedSet<int> foundFullIds;
        Array<int> unmappedFullIds;
        Array<int> activeInputs = getActiveInputs();
        for (int chan : activeInputs) {
            const DataChannel* chanInfo = getDataChannel(chan);
            uint16 sourceNodeId = chanInfo->getSourceNodeID();
            uint16 subProcessorIdx = chanInfo->getSubProcessorIdx();
            int procFullId = int(getProcessorFullId(sourceNodeId, subProcessorIdx));
            foundFullIds.add(procFullId);

            if (subProcessorMap.contains(procFullId)) {
                maxUsedIdx = jmax(maxUsedIdx, subProcessorMap[subProcessorIdx]);
            }
            else {                                                                  // add new entry for this source subprocessor
                if (!subProcessorMap.containsValue(subProcessorIdx)) {  // try to match index if possible
                    subProcessorMap.set(procFullId, subProcessorIdx);
                    maxUsedIdx = jmax(maxUsedIdx, subProcessorIdx);
                }
                else {
                    unmappedFullIds.add(procFullId);
                }
            }
        }

        for (int id : unmappedFullIds) {    // assign remaining unmapped ids
            subProcessorMap.set(id, ++maxUsedIdx);
        }

        Array<int> outdatedFullIds;         // remove outdated entries
        HashMap<int, juce::uint16>::Iterator it(subProcessorMap);
        while (it.next()) {
            int key = it.getKey();
            if (!foundFullIds.contains(key)) {
                outdatedFullIds.add(key);
            }
        }
        for (int id : outdatedFullIds) {
            subProcessorMap.remove(id);
        }
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ASIC VISUALILZATION METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // THSE METHODS ARE USED FOR VISUALIZATION
    void Node::setVisContChan(int newChan) {
        if (newChan >= 0) {
            jassert(newChan < channelInfo.size() && channelInfo[newChan]->isActive());

            // disable event receival temporarily so we can flush the buffer
            int tempVisEventChan = visEventChannel;
            visEventChannel = -1;

            
            // clear timestamp queue
            while (!visASICBuffer->visTsBuffer.empty()) {
                visASICBuffer->visTsBuffer.pop();
            }
            visEventChannel = tempVisEventChan;
        }

        // Clearing the Phase Values
        calVisualizerThreadPtr->clearCompPhaseBuffer();

        // Clearning the Buffer as Well
        visASICBuffer->reset(getDownSample2Rate());

        visContinuousChannel = newChan;

        // If acquisition is stopped (and thus the new channel might be from a different subprocessor),
        // update signal chain. Sinks such as LFP Viewer should receive this information.
        if (!CoreServices::getAcquisitionStatus()) {
            CoreServices::updateSignalChain(getEditor());
        }
    }

    void Node::setVisEventChan(int newChan) {
        jassert(newChan < 8 && newChan >= -1);
        visEventChannel = newChan;
    }


    // miscellaneous helper [Accessed by Canvas.Cpp File]
    int Node::getFullSourceId(int chan)
    {
        const DataChannel* chanInfo = getDataChannel(chan);
        if (!chanInfo)
        {
            jassertfalse;
            return 0;
        }
        uint16 sourceNodeId = chanInfo->getSourceNodeID();
        uint16 subProcessorIdx = chanInfo->getSubProcessorIdx();
        return int(getProcessorFullId(sourceNodeId, subProcessorIdx));
    }


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION THREAD CLASS METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // This Class is Responsible for VISUALIZATION:

    CalVisualizerThread::CalVisualizerThread(Node* nodePtr, int downSample2Rate)
        : Thread("Visualization Thread")
        , node(nodePtr)
        , StopThreadFlag(false)
        , visThreadBuffer(new VisBufferManager(downSample2Rate))
        , visHilbertBufferFFT(new FFTWTransformableArray(int(downSample2Rate * 1.5))) // Equal to Trigger Length {Check: VisBufferManager Class}
        , compPhaseBuffer(new CompPhaseBufferManager()) {
    }

    CalVisualizerThread::~CalVisualizerThread() {
        delete visThreadBuffer;
        delete visHilbertBufferFFT;
        delete compPhaseBuffer;
    }

    void CalVisualizerThread::clearCompPhaseBuffer() {
        compPhaseBuffer->clearBuffer();
    }

    bool CalVisualizerThread::tryToReadVisPhasesfromCompBuffer(std::queue<double>& targetQueue) {
        return compPhaseBuffer->tryToReadVisPhases(targetQueue);
    }


    void CalVisualizerThread::startThisThread() {
        StopThreadFlag = false;    // Setting This Flag Breaks the LOOP of the Thread [We want to get into the loop for Operation]
        startThread();             // Starting the Thread.
        std::cout<<" [VISUALIZER PHASE COMPARISION THREAD!]  The Phase Calculation Thread has been Started"<<std::endl;
    }

    void CalVisualizerThread::stopThisThread() {
        StopThreadFlag = true;         // Set to Break the Thread Loop
        while (isThreadRunning()) {    // If the Thread is Still Running
            waitForThreadToExit(50);   // Wait for 50 milliseconds to exit
        }
        std::cout << " [VISUALIZER PHASE COMPARISION THREAD!]  The Phase Calculation Thread has been Closed Successfully" << std::endl;
    }

    void CalVisualizerThread::run() {
        //std::ofstream threadFileHandle;
        VisBufferManager* visASICBufferPtr = node->getVisBufferPtrASIC();  // Pointer to the main ASIC Buffer node
        double* hilBufferPtr = visHilbertBufferFFT->getRealPointer();      // Getting Real Pointer of Hilbert Buffer


        while (StopThreadFlag == false) {                                    // Run this thread until the end of Acquisition
            if (visASICBufferPtr->isBufferLengthTriggered() && !visASICBufferPtr->visThreadHasCopiedData) {

                visASICBufferPtr->copyBuffer(visThreadBuffer);

                jassert(visHilbertBufferFFT->getLength() == visThreadBuffer->trigIdx); // The Size of FFT buffer should be 1.5 Times DownSampling Rate 2


                // Copying the Contents of the Buffer into the FFTtransformableArray Buffer for Phase Calculation
                Array <float> visHilbertBufferASIC;
                for (int i = 0; i < visASICBufferPtr->trigIdx; i++) {
                    hilBufferPtr[i] = visThreadBuffer->visBPFdata[i];
                    visHilbertBufferASIC.add(visThreadBuffer->visHilbertASIC[i]);
                }

                node->unwrapBuffer(visHilbertBufferASIC.getRawDataPointer(), visThreadBuffer->trigIdx, 0);
                node->smoothBuffer(visHilbertBufferASIC.getRawDataPointer(), visThreadBuffer->trigIdx, 0);

                // Calculating Hilbert Transform;
                visHilbertBufferFFT->hilbert();

                // Will hold the Hilbert FFT Phase Values from HT
                Array <float> mappedAngels;

                // Finding the Angle of [visHilbertBufferFFT] on Each of the event generated by Crossing Detector Plugin
                int visChanIdx           = node->getVisContinuousChannelIndex();       // Visualization Channel Index
                int samplingRate         = node->channelInfo[visChanIdx]->sampleRate;  // 30K
                int downsampledRate      = node->getDownSample2Rate();                 // Downsampling 2 Rate = 1000
                int downsampleRatio      = samplingRate / downsampledRate;             // 30K/1000 = 30
                int LastSampleIndex      = visThreadBuffer->nextIdx;                   // Last Sample Last Idx [Down sampled at 1ms/Sample] 1502
                int trigIndex            = visThreadBuffer->trigIdx;                   // 1500
                int startIndex           = visThreadBuffer->startIdx;                  // 500

                juce::int64 LastSampleTS = visThreadBuffer->lastBuffEndSampTS;         // LastSampleTS [30060]
                juce::int64 minTimeStamp = LastSampleTS - (LastSampleIndex - startIndex) * downsampleRatio;   // The First TS of the Buffer
                juce::int64 maxTimeStamp = LastSampleTS - (LastSampleIndex - trigIndex) * downsampleRatio;

                // Removing the Timestamps That are older than the size [1.5 seconds]  
                while (!visThreadBuffer->visTsBuffer.empty() && visThreadBuffer->visTsBuffer.front() < minTimeStamp) {
                    visThreadBuffer->visTsBuffer.pop();
                }

                // Iterating through the Queue
                while (!visThreadBuffer->visTsBuffer.empty() && visThreadBuffer->visTsBuffer.front() <= maxTimeStamp) {
                    
                    // Getting the TimeStamp
                    juce::int64 ts = visThreadBuffer->visTsBuffer.front();
                    visThreadBuffer->visTsBuffer.pop();
                    
                    // Getting the Hilbert Tansform Phase at the TS
                    int timeStamptoIndx = (ts - minTimeStamp) / downsampleRatio;

                    jassert(timeStamptoIndx >= 0 && timeStamptoIndx <= trigIndex);
                    // Adding startIndx as hilbert Buffer Includes [0-500] Samples
                    float phaseFFT = arg(visHilbertBufferFFT->getAsComplex(timeStamptoIndx + startIndex));
                    mappedAngels.add(phaseFFT);


                    // Passing Event for the Crossing Detector Optimization Plugin [New One]
                    if (!node->getVisPhaseChannel()){
                        jassertfalse; // event channel should not be null here.
                        continue;
                    }

                    double eventData = phaseFFT * 180.0 / Dsp::doublePi; 
                    juce::int64 eventTs = LastSampleTS - node->getNumSamples(visChanIdx);
                    
                    node->addEventThread(LastSampleTS, eventData);

                }

                compPhaseBuffer->insertNewPhaseData(mappedAngels.getRawDataPointer(), mappedAngels.size());


                /*
                // Finding the ANGLES Where the Abrupt Peaks has Occured
                for (int i = visThreadBuffer->startIdx; i < visASICBufferPtr->trigIdx; i++) {

                    if (visHilbertBufferASIC[i - 1] < -165 && visHilbertBufferASIC[i] > 165) {
                        mappedAngels.add(arg(visHilbertBufferFFT->getAsComplex(i))); // Add Radians Angle
                    }
                    else if (visHilbertBufferASIC[i - 1] > 165 && visHilbertBufferASIC[i] < -165) {
                        mappedAngels.add(arg(visHilbertBufferFFT->getAsComplex(i))); // Add Radians Angle
                    }
                }

                compPhaseBuffer->insertNewPhaseData(mappedAngels.getRawDataPointer(), mappedAngels.size());
                */

                //node->editor->addAngleCanvas(mappedAngels);


                //while (!mappedAngels.isEmpty()){
                //    node->editor->canvas

                        //.addAngle(tempPhaseBuffer.front());
                //    mappedAngels.removeLast(1)
                //}


                /*
                // Storing Everything Back into A text File
                threadFileHandle.open("VisualizationTrials.txt", std::ofstream::out | std::ofstream::app);
                for (int i = 0; i < visASICBufferPtr->trigIdx ; i++) {
                    threadFileHandle << i << "," << visHilbertBufferASIC[i] << "," << arg(visHilbertBufferFFT->getAsComplex(i)) * 57.2958 << "," << visThreadBuffer->visBPFdata[i] << std::endl;
                }
                threadFileHandle.close();
                */

            }



        }
    }

	
}


















