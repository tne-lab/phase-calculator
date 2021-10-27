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

#include "PhaseCalculatorCanvas.h"
#include "PhaseCalculatorEditor.h"

namespace PhaseCalculator
{
    Canvas::Canvas(Node* pc)
        : processor(pc)
        , viewport(new Viewport())
        , canvas(new Component("canvas"))
        , rosePlotOptions(new Component("rosePlotOptions"))
        , rosePlot(new RosePlot(this))
    {
        refreshRate = 5;

        canvas->addAndMakeVisible(rosePlot);

        // populate rosePlotOptions
        const Font textFont = Font(18, Font::bold);
        const int textHeight = 25;
        const int indent = 10;
        int xPos = indent;
        int yPos = 10;

        cChannelLabel = new Label("cChannelLabel", "Data channel:");
        cChannelLabel->setBounds(xPos, yPos, 120, textHeight);
        cChannelLabel->setFont(textFont);
        rosePlotOptions->addAndMakeVisible(cChannelLabel);

        eChannelLabel = new Label("eChannelLabel", "Event channel:");
        eChannelLabel->setBounds(xPos += 140, yPos, 120, textHeight);
        eChannelLabel->setFont(textFont);
        rosePlotOptions->addAndMakeVisible(eChannelLabel);

        cChannelBox = new ComboBox("cChannelBox");
        cChannelBox->setTooltip(cChanTooltip);
        cChannelBox->setBounds(xPos = indent + 20, yPos += textHeight + 2, 80, textHeight);
        cChannelBox->addListener(this);
        rosePlotOptions->addAndMakeVisible(cChannelBox);

        eChannelBox = new ComboBox("eChannelBox");
        eChannelBox->setBounds(xPos += 140, yPos, 80, textHeight);
        eChannelBox->addItem("None", 1);
        for (int chan = 1; chan <= 8; ++chan)
        {
            eChannelBox->addItem(String(chan), chan + 1);
        }
        eChannelBox->setSelectedId(1);
        eChannelBox->addListener(this);
        rosePlotOptions->addAndMakeVisible(eChannelBox);

        numBinsLabel = new Label("numBinsLabel", "Number of bins:");
        numBinsLabel->setBounds(xPos = indent, yPos += 2 * textHeight, 200, textHeight);
        numBinsLabel->setFont(textFont);
        rosePlotOptions->addAndMakeVisible(numBinsLabel);

        numBinsSlider = new Slider();
        numBinsSlider->setTextBoxStyle(Slider::TextBoxBelow, false, 40, 30);
        numBinsSlider->setBounds(xPos, yPos += textHeight + 2, optionsWidth - 2 * indent, 40);
        numBinsSlider->setRange(1, RosePlot::maxBins, 1);
        numBinsSlider->setValue(RosePlot::startNumBins);
        numBinsSlider->addListener(this);
        rosePlotOptions->addAndMakeVisible(numBinsSlider);

        clearButton = new UtilityButton("Clear Plot", Font(20, Font::bold));
        clearButton->addListener(this);
        clearButton->setBounds((optionsWidth - 80) / 2, yPos += 55, 80, 30);
        rosePlotOptions->addAndMakeVisible(clearButton);

        referenceLabel = new Label("referenceLabel", "Phase reference:");
        referenceLabel->setBounds(xPos = indent, yPos += 45, 160, textHeight);
        referenceLabel->setFont(textFont);
        rosePlotOptions->addAndMakeVisible(referenceLabel);

        referenceEditable = new Label("referenceEditable", String(RosePlot::startReference));
        referenceEditable->setEditable(true);
        referenceEditable->addListener(rosePlot);
        referenceEditable->setBounds(xPos += 160, yPos, 60, textHeight);
        referenceEditable->setColour(Label::backgroundColourId, Colours::darkgrey);
        referenceEditable->setColour(Label::textColourId, Colours::white);
        referenceEditable->setTooltip(refTooltip);
        rosePlotOptions->addAndMakeVisible(referenceEditable);

        countLabel = new Label("countLabel");
        countLabel->setBounds(xPos = indent, yPos += 45, optionsWidth, textHeight);
        countLabel->setFont(textFont);

        meanLabel = new Label("meanLabel");
        meanLabel->setBounds(xPos, yPos += textHeight, optionsWidth, textHeight);
        meanLabel->setFont(textFont);

        stdLabel = new Label("stdLabel");
        stdLabel->setBounds(xPos, yPos += textHeight, optionsWidth, textHeight);
        stdLabel->setFont(textFont);

        updateStatLabels();
        rosePlotOptions->addAndMakeVisible(countLabel);
        rosePlotOptions->addAndMakeVisible(meanLabel);
        rosePlotOptions->addAndMakeVisible(stdLabel);

        canvas->addAndMakeVisible(rosePlotOptions);

        viewport->setViewedComponent(canvas, false);
        viewport->setScrollBarsShown(true, true);
        addAndMakeVisible(viewport);
    }

    Canvas::~Canvas() {}

    void Canvas::paint(Graphics& g)
    {
        g.fillAll(Colours::grey);
    }

    void Canvas::resized()
    {
        int vpWidth = getWidth();
        int vpHeight = getHeight();
        viewport->setSize(vpWidth, vpHeight);

        int verticalPadding, leftPadding;
        int diameter = getRosePlotDiameter(vpHeight, &verticalPadding);
        int canvasWidth = getContentWidth(vpWidth, diameter, &leftPadding);

        rosePlot->setBounds(leftPadding, verticalPadding, diameter, diameter);

        int optionsX = leftPadding * 2 + diameter;
        int optionsY = verticalPadding + (diameter - minDiameter) / 2;
        rosePlotOptions->setBounds(optionsX, optionsY, optionsWidth, diameter);

        int canvasHeight = diameter + 2 * verticalPadding;
        canvas->setSize(canvasWidth, canvasHeight);
    }

    int Canvas::getRosePlotDiameter(int vpHeight, int* verticalPadding)
    {
        int preferredDiameter = vpHeight - (2 * minPadding);
        int diameter = jmax(minDiameter, jmin(maxDiameter, preferredDiameter));
        if (verticalPadding != nullptr)
        {
            *verticalPadding = jmax(minPadding, (vpHeight - diameter) / 2);
        }
        return diameter;
    }

    int Canvas::getContentWidth(int vpWidth, int diameter, int* leftPadding)
    {
        int widthWithoutPadding = diameter + optionsWidth;
        int lp = jmax(minPadding, jmin(maxLeftPadding, vpWidth - widthWithoutPadding));
        if (leftPadding != nullptr)
        {
            *leftPadding = lp;
        }
        return lp * 2 + widthWithoutPadding;
    }

    void Canvas::refreshState() {}

    void Canvas::update()
    {
        // update continuous channel ComboBox to include only active inputs
        Array<int> activeInputs = processor->getActiveInputs();
        int numActiveInputs = activeInputs.size();
        int numItems = cChannelBox->getNumItems();
        int currSelectedId = cChannelBox->getSelectedId();

        // check whether box needs to be cleared - iterate through existing items and check for correctness
        int startInd = numItems;
        for (int i = 0; i < numItems; ++i)
        {
            if (i >= numActiveInputs ||                           // more items than active inputs
                cChannelBox->getItemId(i) - 1 != activeInputs[i]) // this item doesn't match what it should be
            {
                startInd = 0;
                cChannelBox->clear(dontSendNotification);
            }
        }

        for (int activeChan = startInd; activeChan < numActiveInputs; ++activeChan)
        {
            int chan = activeInputs[activeChan];

            int id = chan + 1;
            cChannelBox->addItem(String(id), id);
            if (id == currSelectedId)
            {
                cChannelBox->setSelectedId(id, dontSendNotification);
            }
        }

        if (cChannelBox->getNumItems() > 0 && cChannelBox->getSelectedId() == 0)
        {
            int firstChannelId = activeInputs[0] + 1;
            cChannelBox->setSelectedId(firstChannelId, dontSendNotification);
        }

        // if channel changed, notify processor of update
        // subtract 1 to change from 1-based to 0-based
        int newId = cChannelBox->getSelectedId();
        if (newId != currSelectedId)
        {
            processor->setParameter(VIS_C_CHAN, static_cast<float>(newId - 1));
        }
    }

    void Canvas::refresh()
    {
        // if no event channel selected, do nothing
        if (eChannelBox->getSelectedId() == 1)
        {
            return;
        }

        // get new angles from visualization phase buffer
        if (processor->tryToReadVisPhases(tempPhaseBuffer))
        {
            // add new angles to rose plot
            while (!tempPhaseBuffer.empty())
            {
                addAngle(tempPhaseBuffer.front());
                tempPhaseBuffer.pop();
            }
        }
    }

    void Canvas::beginAnimation()
    {
        // disable continuous channel options from a different subprocessor
        // this is necessary to maintain the correct source subprocessor metadata on the
        // visualized phase event channel, to avoid problems with the LFP viewer.
        // (shouldn't cause issues unless the setup is pretty weird)
        int selectedId = cChannelBox->getSelectedId();
        if (selectedId > 0)
        {
            int currSourceSubproc = processor->getFullSourceId(selectedId - 1);

            int numItems = cChannelBox->getNumItems();
            for (int i = 0; i < numItems; ++i)
            {
                int id = cChannelBox->getItemId(i);
                if (id != selectedId && processor->getFullSourceId(id - 1) != currSourceSubproc)
                {
                    cChannelBox->setItemEnabled(id, false);
                }
            }
        }

        startCallbacks();
    }

    void Canvas::endAnimation()
    {
        // re-enable all continuous channel options
        int numItems = cChannelBox->getNumItems();
        for (int i = 0; i < numItems; ++i)
        {
            cChannelBox->setItemEnabled(cChannelBox->getItemId(i), true);
        }

        stopCallbacks();
    }

    void Canvas::setParameter(int, float) {}
    void Canvas::setParameter(int, int, int, float) {}

    void Canvas::addAngle(double newAngle)
    {
        rosePlot->addAngle(newAngle);
        updateStatLabels();
    }

    void Canvas::clearAngles()
    {
        rosePlot->clear();
        updateStatLabels();
    }

    void Canvas::comboBoxChanged(ComboBox* comboBoxThatHasChanged)
    {
        if (comboBoxThatHasChanged == cChannelBox)
        {
            // subtract 1 to change from 1-based to 0-based
            float newValue = static_cast<float>(cChannelBox->getSelectedId() - 1);
            processor->setParameter(VIS_C_CHAN, newValue);
        }
        else if (comboBoxThatHasChanged == eChannelBox)
        {
            // subtract 2, since index 1 == no channel (-1)
            float newValue = static_cast<float>(eChannelBox->getSelectedId() - 2);
            processor->setParameter(VIS_E_CHAN, newValue);
        }
    }

    void Canvas::sliderValueChanged(Slider* slider)
    {
        if (slider == numBinsSlider)
        {
            rosePlot->setNumBins(static_cast<int>(slider->getValue()));
        }
    }

    void Canvas::buttonClicked(Button* button)
    {
        if (button == clearButton)
        {
            clearAngles();
        }
    }

    void Canvas::updateStatLabels()
    {
        int numAngles = rosePlot->getNumAngles();
        double mean = std::round(100 * rosePlot->getCircMean()) / 100;
        double stddev = std::round(100 * rosePlot->getCircStd()) / 100;

        countLabel->setText(String::formatted(countFmt, numAngles), dontSendNotification);
        meanLabel->setText(String::formatted(meanFmt, mean), dontSendNotification);
        stdLabel->setText(String::formatted(stdFmt, stddev), dontSendNotification);
    }

    void Canvas::saveVisualizerParameters(XmlElement* xml)
    {
        XmlElement* visValues = xml->createNewChildElement("VISUALIZER");
        visValues->setAttribute("eventChannelId", eChannelBox->getSelectedId());
        visValues->setAttribute("numBins", numBinsSlider->getValue());
        visValues->setAttribute("phaseRef", referenceEditable->getText());
    }

    void Canvas::loadVisualizerParameters(XmlElement* xml)
    {
        forEachXmlChildElementWithTagName(*xml, xmlNode, "VISUALIZER")
        {
            int eventChannelId = xmlNode->getIntAttribute("eventChannelId", eChannelBox->getSelectedId());
            if (eChannelBox->indexOfItemId(eventChannelId) != -1)
            {
                eChannelBox->setSelectedId(eventChannelId, sendNotificationSync);
            }

            numBinsSlider->setValue(xmlNode->getDoubleAttribute("numBins", numBinsSlider->getValue()), sendNotificationSync);
            referenceEditable->setText(xmlNode->getStringAttribute("phaseRef", referenceEditable->getText()), sendNotificationSync);
        }
    }

    void Canvas::displayContinuousChan(int chan)
    {
        if (cChannelBox->indexOfItemId(chan + 1) == -1)
        {
            jassertfalse;
            return;
        }
        // remember to switch to 1-based
        cChannelBox->setSelectedId(chan + 1, dontSendNotification);
    }

    /**** RosePlot ****/

    RosePlot::RosePlot(Canvas* c)
        : canvas(c)
        , referenceAngle(static_cast<double>(startReference))
        , numBins(startNumBins)
        , faceColor(Colour(230, 168, 0))
        , edgeColor(Colours::black)
        , bgColor(Colours::black)
        , edgeWeight(1)
        , rSum(0)
    {
        updateAngles();
        reorganizeAngleData();
    }

    RosePlot::~RosePlot() {}

    void RosePlot::paint(Graphics& g)
    {
        // dimensions
        juce::Rectangle<float> bounds = getBounds().toFloat();
        float squareSide = jmin(bounds.getHeight(), bounds.getWidth() - 2 * textBoxSize);
        juce::Rectangle<float> plotBounds = bounds.withZeroOrigin();
        plotBounds = plotBounds.withSizeKeepingCentre(squareSide, squareSide);
        g.setColour(bgColor);
        g.fillEllipse(plotBounds);

        // draw grid
        // spokes and degree labels (every 30 degrees)
        Point<float> center = plotBounds.getCentre();
        Line<float> spoke(center, center);
        juce::Rectangle<int> textBox(textBoxSize, textBoxSize);
        g.setFont(Font(textBoxSize / 2, Font::bold));
        for (int i = 0; i < 12; ++i)
        {
            float juceAngle = i * float_Pi / 6;
            spoke.setEnd(center.getPointOnCircumference(squareSide / 2, juceAngle));
            g.setColour(Colours::lightgrey);
            g.drawLine(spoke);

            float textRadius = (squareSide + textBoxSize) / 2;
            Point<int> textCenter = center.getPointOnCircumference(textRadius, juceAngle).toInt();
            int degreeAngle = (450 - 30 * i) % 360;
            g.setColour(Colours::black);
            g.drawFittedText(String(degreeAngle), textBox.withCentre(textCenter), Justification::centred, 1);
        }

        // concentric circles
        int nCircles = 3;
        juce::Rectangle<float> circleBounds;
        g.setColour(Colours::lightgrey);
        for (int i = 1; i < nCircles; ++i)
        {
            float diameter = (squareSide * i) / nCircles;
            circleBounds = plotBounds.withSizeKeepingCentre(diameter, diameter);
            g.drawEllipse(circleBounds, 1);
        }

        // get count for each rose plot segment
        int nSegs = binMidpoints.size();
        Array<int> segmentCounts;
        segmentCounts.resize(nSegs);
        int maxCount = 0;
        int totalCount = 0;
        for (int seg = 0; seg < nSegs; ++seg)
        {
            int count = static_cast<int>(angleData->count(binMidpoints[seg] + referenceAngle));
            segmentCounts.set(seg, count);
            maxCount = jmax(maxCount, count);
            totalCount += count;
        }

        jassert(totalCount == angleData->size());
        jassert((maxCount == 0) == (angleData->empty()));

        // construct path
        Path rosePath;
        for (int seg = 0; seg < nSegs; ++seg)
        {
            if (segmentCounts[seg] == 0)
            {
                continue;
            }

            float size = squareSide * segmentCounts[seg] / static_cast<float>(maxCount);
            rosePath.addPieSegment(plotBounds.withSizeKeepingCentre(size, size),
                segmentAngles[seg].first, segmentAngles[seg].second, 0);
        }

        // paint path
        g.setColour(faceColor);
        g.fillPath(rosePath);
        g.setColour(edgeColor);
        g.strokePath(rosePath, PathStrokeType(edgeWeight));
    }

    void RosePlot::setNumBins(int newNumBins)
    {
        if (newNumBins != numBins && newNumBins > 0 && newNumBins <= maxBins)
        {
            numBins = newNumBins;
            updateAngles();
            reorganizeAngleData();
            repaint();
        }
    }

    void RosePlot::setReference(double newReference)
    {
        if (newReference != referenceAngle)
        {
            referenceAngle = newReference;
            reorganizeAngleData();
            repaint();
        }
    }

    void RosePlot::addAngle(double newAngle)
    {
        newAngle = Node::circDist(newAngle, 0.0);
        angleData->insert(newAngle);
        rSum += std::exp(std::complex<double>(0, newAngle));
        repaint();
    }

    void RosePlot::clear()
    {
        angleData->clear();
        rSum = 0;
        repaint();
    }

    int RosePlot::getNumAngles()
    {
        return static_cast<int>(angleData->size());
    }

    double RosePlot::getCircMean(bool usingReference)
    {
        if (angleData->empty())
        {
            return 0;
        }

        double reference = usingReference ? referenceAngle : 0.0;
        // use range of (-90, 270] for ease of use
        double meanRad = Node::circDist(std::arg(rSum), reference, 3 * double_Pi / 2);
        return radiansToDegrees(meanRad);
    }

    double RosePlot::getCircStd()
    {
        if (angleData->empty())
        {
            return 0;
        }

        double r = std::abs(rSum) / angleData->size();
        double stdRad = std::sqrt(-2 * std::log(r));
        return radiansToDegrees(stdRad);
    }

    void RosePlot::labelTextChanged(Label* labelThatHasChanged)
    {
        if (labelThatHasChanged->getName() == "referenceEditable")
        {
            double doubleInput;
            double currReferenceDeg = radiansToDegrees(referenceAngle);
            bool valid = Editor::updateControl(labelThatHasChanged,
                -DBL_MAX, DBL_MAX, currReferenceDeg, doubleInput);

            if (valid)
            {
                // convert to radians
                double newReference = Node::circDist(degreesToRadians(doubleInput), 0.0);
                labelThatHasChanged->setText(String(radiansToDegrees(newReference)), dontSendNotification);
                setReference(newReference);
                canvas->updateStatLabels();
            }
        }
    }


    /*** RosePlot private members ***/

    RosePlot::AngleDataMultiset::AngleDataMultiset(int numBins, double referenceAngle)
        : std::multiset<double, std::function<bool(double, double)>>(
        std::bind(circularBinCompare, numBins, referenceAngle, std::placeholders::_1, std::placeholders::_2))
    {}

    RosePlot::AngleDataMultiset::AngleDataMultiset(int numBins, double referenceAngle, AngleDataMultiset* dataSource)
        : AngleDataMultiset(numBins, referenceAngle)
    {
        insert(dataSource->begin(), dataSource->end());
    }


    bool RosePlot::AngleDataMultiset::circularBinCompare(int numBins, double referenceAngle, double lhs, double rhs)
    {
        double lhsDist = Node::circDist(lhs, referenceAngle);
        double rhsDist = Node::circDist(rhs, referenceAngle);
        int lhsBin = int(lhsDist * numBins / (2 * double_Pi));
        int rhsBin = int(rhsDist * numBins / (2 * double_Pi));
        return lhsBin < rhsBin;
    }

    void RosePlot::reorganizeAngleData()
    {
        ScopedPointer<AngleDataMultiset> newAngleData;
        if (angleData == nullptr)
        {
            // construct empty container
            newAngleData = new AngleDataMultiset(numBins, referenceAngle);
        }
        else
        {
            // copy existing data to new container
            newAngleData = new AngleDataMultiset(numBins, referenceAngle, angleData.get());
        }

        angleData.swapWith(newAngleData);
    }

    void RosePlot::updateAngles()
    {
        float step = 2 * float_Pi / numBins;
        binMidpoints.resize(numBins);
        for (int i = 0; i < numBins; ++i)
        {
            binMidpoints.set(i, step * (i + 0.5));
            float firstAngle = float(Node::circDist(double_Pi / 2, step * (i + 1)));
            segmentAngles.set(i, { firstAngle, firstAngle + step });
        }
    }
}