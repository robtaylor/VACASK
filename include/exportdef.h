#ifndef __EXPORTDEF_DEFINED
#define __EXPORTDEF_DEFINED

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_SIM_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #elseif USING_SIM_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #else
    #define DLL_PUBLIC
  #endif
  #define DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #if defined BUILDING_SIM_DLL || defined USING_SIM_DLL
      #define DLL_PUBLIC __attribute__ ((visibility ("default")))
      #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
    #else
      #define DLL_PUBLIC
      #define DLL_LOCAL
    #endif
  #else
    #define DLL_PUBLIC
    #define DLL_LOCAL
  #endif
#endif

#endif
