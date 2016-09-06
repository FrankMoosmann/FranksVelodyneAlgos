/*!
    \file   Deprecated
    \brief  Provides a typedef to declare functions deprecated
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    1.2.2010

    Copyright: Frank Moosmann
*/

#ifndef DEPRECATED_H_
#define DEPRECATED_H_

// Each compiler has its own deprecated-warning
#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif
// You will encounter problems if a function return type has a commas in its name
// e.g. std::pair<int, int> as this will be interpreted by the preprocesor as passing
// 2 arguments to the DEPRECATED macro.
// In that case you would have to typedef the return type.

// example usage:
// DEPRECATED(void OldFunc(int a, float b)); // don't use me any more
// void NewFunc(int a, double b); // use me instead


#endif /* DEPRECATED_H_ */
