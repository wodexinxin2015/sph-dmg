/************************************************************************************
SPH Discrete Model Generator
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#pragma once

#ifdef __ANDROID__
#define win32 0
#elif __linux__
#define win32 0
#elif _WIN32
#define win32 1
#endif

#define size 5000000