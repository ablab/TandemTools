/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/

#ifndef __JELLYFISH_ATOMIC_GCC_HPP__
#define __JELLYFISH_ATOMIC_GCC_HPP__

namespace atomic
{
  class gcc
  {
  public:
    template<typename T>
    static inline T cas(volatile T *ptr, T oval, T nval) {
      return __sync_val_compare_and_swap(ptr, oval, nval);
    }

    template<typename T>
    static inline T set(T *ptr, T nval) {
      return __sync_lock_test_and_set(ptr, nval);
    }

    template<typename T>
    static inline T add_fetch(volatile T *ptr, T x) {
      T ncount = *ptr, count;
      do {
	count = ncount;
	ncount = cas((T *)ptr, count, count + x);
      } while(ncount != count);
      return count + x;
    }

    template<typename T>
    static inline T fetch_add(volatile T *ptr, T x) {
      T ncount = *ptr, count;
      do {
	count = ncount;
	ncount = cas((T *)ptr, count, (T)(count + x));
      } while(ncount != count);
      return count;
    }

    template<typename T>
    static inline T set_to_max(volatile T *ptr, T x) {
      T count = *ptr;
      while(x > count) {
        T ncount = cas(ptr, count, x);
        if(ncount == count)
          return x;
        count = ncount;
      }
      return count;
    }
  };
}
#endif
