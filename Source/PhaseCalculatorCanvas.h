/*
------------------------------------------------------------------

This file is part of a plugin for the Open Ephys GUI
Copyright (C) 2018 Translational NeuroEngineering Laboratory, MGH

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

#ifndef PHASE_CALCULATOR_CANVAS_H_INCLUDED
#define PHASE_CALCULATOR_CANVAS_H_INCLUDED

#include "PhaseCalculator.h"
#include <VisualizerWindowHeaders.h>
#include <set> // std::multiset

class PhaseCalculatorCanvas;

class RosePlot
    : public Component
    , public Label::Listener
{
public:
    RosePlot(PhaseCalculatorCanvas* c);
    ~RosePlot();

    void paint(Graphics& g) override;
    
    // Change number of bins and repaint
    void setNumBins(int newNumBins);

    // Change reference angle and repaint
    void setReference(double newReference);

    // Add a new angle (in radians) and repaint
    void addAngle(double newAngle);

    // Remove all angles from the plot and repaint
    void clear();

    int getNumAngles();

    // output statistics, in degrees
    double getCircMean(bool usingReference = true);
    double getCircStd();

    // implements Label::Listener
    void labelTextChanged(Label* labelThatHasChanged) override;

    static const int MAX_BINS = 120;
    static const int START_NUM_BINS = 24;
    static const int START_REFERENCE = 0;
    static const int TEXT_BOX_SIZE = 50;

private:
    struct circularBinComparator
    {
        circularBinComparator(int numBinsIn, double referenceAngleIn);

        // Depending on numBins and reference, returns true iff lhs belongs in a lower bin than rhs.
        bool operator() (const double& lhs, const double& rhs) const;

    private:
        int numBins;
        double referenceAngle;
    };

    class AngleDataMultiset : public std::multiset<double, circularBinComparator>
    {
    public:
        // create empty multiset
        AngleDataMultiset(int numBins, double referenceAngle);

        // copy nodes from dataSource to newly constructed multiset
        AngleDataMultiset(int numBins, double referenceAngle, AngleDataMultiset* dataSource);
    };

    // make binMidpoints and segmentAngles reflect current numBins
    void updateAngles();

    // reassign angleData to refer to a new AngleDataMultiset with the current parameters
    // (but keep the same data points)
    void reorganizeAngleData();

    PhaseCalculatorCanvas* canvas;

    ScopedPointer<AngleDataMultiset> angleData;
    int numBins;
    double referenceAngle;

    // for each rose plot segment:
    // midpoint angle in radians CCW from positive x axis
    Array<double> binMidpoints;
    // inputs to addPieSegment (clockwise from top)
    Array<std::pair<float, float>> segmentAngles;

    Colour faceColor;
    Colour edgeColor;
    Colour bgColor;
    float edgeWeight;

    // sum of exp(a*j) for each angle a
    std::complex<double> rSum;

    static const double PI;
    static const float PI_F;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RosePlot);
};

class PhaseCalculatorCanvas 
    : public Visualizer
    , public ComboBox::Listener
    , public Slider::Listener
    , public Button::Listener
{
public:
    PhaseCalculatorCanvas(PhaseCalculator* pc);
    ~PhaseCalculatorCanvas();
    void refreshState() override;
    void update() override;
    void refresh() override;
    void beginAnimation() override;
    void endAnimation() override;
    void setParameter(int, float) override;
    void setParameter(int, int, int, float) override;

    void paint(Graphics& g) override;
    void resized() override;

    void addAngle(double newAngle);
    void clearAngles();

    // implements ComboBox::Listener
    void comboBoxChanged(ComboBox* comboBoxThatHasChanged) override;

    // implements Slider::Listener
    void sliderValueChanged(Slider* slider) override;

    // implements Button::Listener
    void buttonClicked(Button* button) override;

    // updates countLabel, meanLabel, and stdLabel
    void updateStatLabels();

    void saveVisualizerParameters(XmlElement* xml) override;
    void loadVisualizerParameters(XmlElement* xml) override;

    // display updaters - do not trigger listeners.
    void displayContinuousChan(int chan);

private:
    int getRosePlotDiameter(int height, int* verticalPadding);
    int getContentWidth(int width, int diameter, int* leftPadding);

    PhaseCalculator* processor;

    ScopedPointer<Viewport>  viewport;
    ScopedPointer<Component> canvas;
    ScopedPointer<Component> rosePlotOptions;
    ScopedPointer<RosePlot>  rosePlot;

    // options panel
    ScopedPointer<Label>     cChannelLabel;
    ScopedPointer<ComboBox>  cChannelBox;
    ScopedPointer<Label>     eChannelLabel;
    ScopedPointer<ComboBox>  eChannelBox;

    ScopedPointer<Label>  numBinsLabel;
    ScopedPointer<Slider> numBinsSlider;

    ScopedPointer<UtilityButton> clearButton;

    ScopedPointer<Label> referenceLabel;
    ScopedPointer<Label> referenceEditable;
    
    ScopedPointer<Label> countLabel;
    ScopedPointer<Label> meanLabel;
    ScopedPointer<Label> stdLabel;

    static const int MIN_PADDING = 5;
    static const int MAX_LEFT_PADDING = 50;
    static const int MIN_DIAMETER = 350;
    static const int MAX_DIAMETER = 550;
    static const int OPTIONS_WIDTH = 320;

    const String C_CHAN_TOOLTIP = "Channel containing data whose high-accuracy phase is calculated for each event";
    const String REF_TOOLTIP = "Base phase (in degrees) to subtract from each calculated phase";
    const String COUNT_FMT = L"Events received: %d";
    const String MEAN_FMT = L"Mean phase (vs. reference): %.2f\u00b0";
    const String STD_FMT = L"Standard deviation phase: %.2f\u00b0";

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseCalculatorCanvas);
};

#endif // PHASE_CALCULATOR_CANVAS_H_INCLUDED