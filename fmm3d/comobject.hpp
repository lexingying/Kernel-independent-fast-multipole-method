/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#ifndef _COMOBJECT_HPP_
#define _COMOBJECT_HPP_

#include "common/commoninc.hpp"

using std::string;
using std::map;

class ComObject
{
protected:
  string _prefix;
public:
  ComObject(const string& prefix): _prefix(prefix) {;}
  virtual ~ComObject() {;}
  //-------------------------
  const string& prefix() { return _prefix; }  //virtual int setFromOptions(map<string,string>&) { return 0; }
};

#endif
