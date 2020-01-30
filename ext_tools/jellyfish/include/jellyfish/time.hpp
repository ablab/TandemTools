/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/

#ifndef __JELLYFISH_TIME_HPP__
#define __JELLYFISH_TIME_HPP__

#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <iomanip>

class Time {
  static const suseconds_t max_useconds = 1000000UL;
  struct timeval tv;

 public:
  static const Time zero;
  explicit Time(bool init = true) {
    if(init)
      now();
  }
  Time(time_t sec, suseconds_t usec) {
    tv.tv_sec  = sec;
    tv.tv_usec = usec;
  }
  Time &operator=(const Time &o) {
    if(&o != this) {
      tv.tv_sec = o.tv.tv_sec;
      tv.tv_usec = o.tv.tv_usec;
    }
    return *this;
  }

  Time & operator-=(const Time &o) {
    tv.tv_sec -= o.tv.tv_sec;
    if(o.tv.tv_usec > tv.tv_usec) {
      tv.tv_usec = (max_useconds + tv.tv_usec) - o.tv.tv_usec;
      --tv.tv_sec;
    } else {
      tv.tv_usec -= o.tv.tv_usec;
    }
    return *this;
  }
  const Time operator-(const Time &o) const {
    return Time(*this) -= o;
  }

  Time & operator+=(const Time &o) {
    tv.tv_sec  += o.tv.tv_sec;
    tv.tv_usec += o.tv.tv_usec;
    if(tv.tv_usec >= max_useconds) {
      ++tv.tv_sec;
      tv.tv_usec -= max_useconds;
    }
    return *this;
  }
  const Time operator+(const Time &o) const {
    return Time(*this) += o;
  }

  bool operator<(const Time& o) const {
    return tv.tv_sec < o.tv.tv_sec || (tv.tv_sec == o.tv.tv_sec && tv.tv_usec < o.tv.tv_usec);
  }

  void now() { gettimeofday(&tv, NULL); }
  Time elapsed() const {
    return Time() - *this;
  }


  std::string str() const {
    std::ostringstream res;
    res << tv.tv_sec << "."
        << std::setfill('0') << std::setw(6) << std::right << tv.tv_usec;
    return res.str();
  }
};

#endif // __TIME_HPP__
