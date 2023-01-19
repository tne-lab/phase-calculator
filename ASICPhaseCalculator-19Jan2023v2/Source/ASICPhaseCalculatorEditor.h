/*
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


#ifndef ASIC_PHASE_CALCULATOR_EDITOR_H_INCLUDED
#define ASIC_PHASE_CALCULATOR_EDITOR_H_INCLUDED

#include <VisualizerEditorHeaders.h>
#include "ASICPhaseCalculator.h"
#include "ASICPhaseCalculatorCanvas.h"


namespace ASICPhaseCalculator {


	class Editor : public VisualizerEditor {
		
	public:

		friend class AllChangeListener;

		/** Constructor */
		Editor(Node* parentNode);

		/** Destructor */
		~Editor();

		// update display based on current channel
		void channelChanged(int chan, bool newState) override;
		
		std::unique_ptr<Node> node;                   // A Pointer Declaration which will point to the ASIC Processor  
		std::unique_ptr<AllChangeListener> chglistnr; // A Pointer Declaration to the allChangeListener Class
		std::unique_ptr<Canvas> canvas;               // A Pointer Declaration to the Chanvas Object

		// enable/disable controls when acquisiton starts/ends
		void startAcquisition() override;
		void stopAcquisition() override;
		
		/** Generating Canvas*/
		Visualizer* createNewCanvas() override;

	private:


		/** Return whether the control contained a valid input between min and max, inclusive. If so, it is stored 
		    in *out and the control is updated with the parsed input. Otherwise, the control is reset to defaultValue.
            In header to make sure specializations not used in ASICPhaseCalculatorEditor.cpp
            are still available to other translation units. */
		template<typename Ctrl, typename T>
		static bool updateControl(Ctrl* c, const T min, const T max, const T defaultValue, T& out) {
			if (readNumber(c->getText(), out)) {
				out = jmax(min, jmin(max, out));
				c->setText(String(out), dontSendNotification);
				return true;
			}
			c->setText(String(defaultValue), dontSendNotification);
			return false;
		}
		

		/** Tries to read a number of type out from input. Returns false if unsuccessful.
		    Otherwise, returns true and writes the result to *out. */
		template<typename T>
		static bool readNumber(const String& input, T& out){
			std::istringstream istream(input.toStdString());
			istream >> out;
			return !istream.fail();
		}


		ScopedPointer<Label>    bandLabel;
		ScopedPointer<ComboBox> bandBox;
		ScopedPointer<Label>    outputModeLabel;
		ScopedPointer<ComboBox> outputModeBox;
		
		//added by sumedh
		ScopedPointer<Label> downsampleIlabel;
		ScopedPointer<Label> downsampleIIlabel;
		ScopedPointer<Label> savedDataTxtlabel;
		ScopedPointer<ComboBox> istextsavebox;
		ScopedPointer<ComboBox> isLAAbox;
		ScopedPointer<Label>    downsampleIlabelEdit;
		ScopedPointer<Label>    downsampleIIlabelEdit;
		ScopedPointer<Label>    txtlabelTEXT;
		ScopedPointer<Label>    txtlabelLAA;


		// added by Sumedh
		// constants
		const String freqRangeTooltip = "Each option corresponds internally to a Hilbert transformer " +
			String("that is optimized for this frequency range. After selecting a range, you can adjust ") +
			"'low' and 'high' to filter to any passband within this range.";
		const String outputModeTooltip = "Which component of the analytic signal to output. If 'PH+MAG' is selected, " +
			String("creates a second channel for each enabled input channel and outputs phases ") +
			"on the original channels and magnitudes on the corresponding new channels.";

		// added by Sumedh
		const String TxtTooltip = "Select YES to save Data into a text file";
		const String LAATooltip = "Select YES to enable LAA";

		/** Generates an assertion if this class leaks */
		JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Editor);
	};
	



	 
	

	// allChangeListener is the DropDown Update Listener Class
	class AllChangeListener
		: public ComboBox::Listener
		, public Label::Listener {
	public:

		AllChangeListener(Editor* editorNode);
		~AllChangeListener();

		/** This is Combo-Box Listener Method, Whever "comboBox" is updated, This method Execultes
			We override this listener Menthod for doing the Necessary Operations When executed */
		void comboBoxChanged(ComboBox* comboBoxThatHasChanged) override;

		/** This is Edit Text Listener Method, Whever the Edit TxT box is updated, This method Executes
		    We override this method to do the necessary Operations based on Our Requirement*/
		void labelTextChanged(Label* labelThatHasChanged) override; 

	private:
		std::unique_ptr<Editor> editor;  // Pointer to the ASIC Editor Class
		JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AllChangeListener);
	};


}



#endif // PROCESSORPLUGINEDITOR_H_DEFINED