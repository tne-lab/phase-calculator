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


#include <PluginInfo.h>

#include "ASICPhaseCalculator.h"
#include <string>

#ifdef WIN32
#include <Windows.h>
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

using namespace Plugin;

#define NUM_PLUGINS 1

extern "C" EXPORT void getLibInfo(Plugin::LibraryInfo* info)
{
	/* API version, defined by the GUI source.
	Should not be changed to ensure it is always equal to the one used in the latest codebase.
	The GUI refueses to load plugins with mismatched API versions */
	info->apiVersion = PLUGIN_API_VER;
	info->name = "ASIC Phase Calculator"; // <---- update
	info->libVersion = 1; // <---- update
	info->numPlugins = NUM_PLUGINS;
}

extern "C" EXPORT int getPluginInfo(int index, Plugin::PluginInfo* info)
{
	switch (index)
	{
		//one case per plugin. This example is for a processor which connects directly to the signal chain
	case 0:
		//Type of plugin. See "Source/Processors/PluginManager/OpenEphysPlugin.h" for complete info about the different type structures
		info->type = PLUGIN_TYPE_PROCESSOR;

		//Processor name
		info->processor.name = "ASIC"; //Processor name shown in the GUI

		//Type of processor. Can be FILTER, SOURCE, SINK or UTILITY. Specifies where on the processor list will appear
		info->processor.type = Plugin::FilterProcessor;

		//Class factory pointer. Replace "ProcessorPluginSpace::ProcessorPlugin" with the namespace and class name.
		info->processor.creator = &(Plugin::createProcessor<ASICPhaseCalculator::Node>);
		break;
	default:
		return -1;
		break;
	}
	return 0;
}

#ifdef WIN32
BOOL WINAPI DllMain(IN HINSTANCE hDllHandle,
	IN DWORD     nReason,
	IN LPVOID    Reserved)
{
	return TRUE;
}

#endif
