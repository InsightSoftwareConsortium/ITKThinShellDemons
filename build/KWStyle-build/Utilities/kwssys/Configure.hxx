/* Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
   file Copyright.txt or https://cmake.org/licensing#kwsys for details.  */
#ifndef kwssys_Configure_hxx
#define kwssys_Configure_hxx

/* Include C configuration.  */
#include <kwssys/Configure.h>

/* Whether wstring is available.  */
#define kwssys_STL_HAS_WSTRING 1
/* Whether <ext/stdio_filebuf.h> is available. */
#define kwssys_CXX_HAS_EXT_STDIO_FILEBUF_H                         \
  1
/* Whether the translation map is available or not. */
#define kwssys_SYSTEMTOOLS_USE_TRANSLATION_MAP                     \
  1

#if defined(__SUNPRO_CC) && __SUNPRO_CC > 0x5130 && defined(__has_attribute)
#  define kwssys__has_cpp_attribute(x) __has_attribute(x)
#elif defined(__has_cpp_attribute)
#  define kwssys__has_cpp_attribute(x) __has_cpp_attribute(x)
#else
#  define kwssys__has_cpp_attribute(x) 0
#endif

#if __cplusplus >= 201103L
#  define kwssys_NULLPTR nullptr
#else
#  define kwssys_NULLPTR 0
#endif

#ifndef kwssys_FALLTHROUGH
#  if __cplusplus >= 201703L &&                                               \
    kwssys__has_cpp_attribute(fallthrough)
#    define kwssys_FALLTHROUGH [[fallthrough]]
#  elif __cplusplus >= 201103L &&                                             \
    kwssys__has_cpp_attribute(gnu::fallthrough)
#    define kwssys_FALLTHROUGH [[gnu::fallthrough]]
#  elif __cplusplus >= 201103L &&                                             \
    kwssys__has_cpp_attribute(clang::fallthrough)
#    define kwssys_FALLTHROUGH [[clang::fallthrough]]
#  endif
#endif
#ifndef kwssys_FALLTHROUGH
#  define kwssys_FALLTHROUGH static_cast<void>(0)
#endif

#undef kwssys__has_cpp_attribute

/* If building a C++ file in kwsys itself, give the source file
   access to the macros without a configured namespace.  */
#if defined(KWSYS_NAMESPACE)
#  if !kwssys_NAME_IS_KWSYS
#    define kwsys kwssys
#  endif
#  define KWSYS_NAME_IS_KWSYS kwssys_NAME_IS_KWSYS
#  define KWSYS_STL_HAS_WSTRING kwssys_STL_HAS_WSTRING
#  define KWSYS_CXX_HAS_EXT_STDIO_FILEBUF_H                                   \
    kwssys_CXX_HAS_EXT_STDIO_FILEBUF_H
#  define KWSYS_FALLTHROUGH kwssys_FALLTHROUGH
#  define KWSYS_SYSTEMTOOLS_USE_TRANSLATION_MAP                               \
    kwssys_SYSTEMTOOLS_USE_TRANSLATION_MAP
#endif

#endif
