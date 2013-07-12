////////////////////////////////////////////////////////////////////////////////
//  
//  Copyright (C) 2012 Eric Conte, Benjamin Fuks, Guillaume Serret
//  The MadAnalysis development team, email: <ma5team@iphc.cnrs.fr>
//  
//  This file is part of MadAnalysis 5.
//  Official website: <http://madanalysis.irmp.ucl.ac.be>
//  
//  MadAnalysis 5 is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  
//  MadAnalysis 5 is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with MadAnalysis 5. If not, see <http://www.gnu.org/licenses/>
//  
////////////////////////////////////////////////////////////////////////////////


#ifndef GZ_ISTREAM_H
#define GZ_ISTREAM_H

// STL headers
#include <iostream>
#include <fstream>

// ZLib header
#include <zlib.h>


// -------------------------------------------------------------
//                   CLASS IGZ_ISTREAMBUF
// -------------------------------------------------------------
class gz_istreambuf : public std::streambuf
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:

  /// size of the data buffer
  static const int bufferSize = 47+256;    

  gzFile file;               // file handle for compressed file
  char   buffer[bufferSize]; // data buffer
  char   opened;             // open/close state of stream
  int    mode;               // I/O mode


  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 private :

  /// Flush the buffer
  int flush_buffer();

 public:

  /// Constructor withtout arguments
  gz_istreambuf() : opened(0)
  {
    setp( buffer, buffer + (bufferSize-1));
    setg( buffer + 4, buffer + 4, buffer + 4); 
  }

  /// Is opened
  int is_open() { return opened; }

  /// Opening the gzip file
  gz_istreambuf* open( const char* name, int open_mode);

  /// Closing the file
  gz_istreambuf* close();

  /// Destructor
  ~gz_istreambuf() { close(); }
    
  /// Overflow
  virtual int overflow( int c = EOF);

  /// Underflow
  virtual int underflow();

  /// Synchronize input buffer
  virtual int sync();
};


// -------------------------------------------------------------
//                   CLASS IGZ_ISTREAMBASE
// -------------------------------------------------------------
class gz_istreambase : virtual public std::ios
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 protected:
  gz_istreambuf buf;


  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 public:

  /// Constructor without arguments
  gz_istreambase() 
  { init(&buf); }

  /// Constructor with arguments
  gz_istreambase( const char* name, int open_mode)
  {
    init( &buf);
    open( name, open_mode);
  }

  /// Destructor
  ~gz_istreambase()
  { buf.close(); }

  /// Open a gzip file
  void open( const char* name, int open_mode)
  {
    if (!buf.open( name, open_mode))
        clear( rdstate() | std::ios::badbit);
  }

  /// Close a gzip file
  void close()
  {
    if (buf.is_open())
    {
      if (!buf.close())
        clear( rdstate() | std::ios::badbit);
    }
  }

  /// Read the buffer
  gz_istreambuf* rdbuf()
  { return &buf; }

};



// -------------------------------------------------------------
//                      CLASS IGZ_ISTREAM
// -------------------------------------------------------------

class gz_istream : public gz_istreambase, public std::istream
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------

  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 public:


  /// Constructor without arguments
  gz_istream() : std::istream(&buf)
  {}

  /// Constructor with arguments
  gz_istream( const char* name, int open_mode = std::ios::in)
      : gz_istreambase( name, open_mode), std::istream( &buf) 
  {}
  
  /// Read buffer
  gz_istreambuf* rdbuf()
  { return gz_istreambase::rdbuf(); }

  /// open a gzip file
  void open( const char* name, int open_mode = std::ios::in)
  { gz_istreambase::open( name, open_mode); }

};


#endif

