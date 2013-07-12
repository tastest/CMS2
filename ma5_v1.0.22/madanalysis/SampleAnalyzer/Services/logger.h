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


#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <sstream>


#define DEBUG   Logger::getInstance()->debug()
#define INFO    Logger::getInstance()->info()
#define WARNING Logger::getInstance()->warning(__FILE__,__LINE__)
#define ERROR   Logger::getInstance()->error(__FILE__,__LINE__)


class LoggerStream
{
public:

  enum ColorType {NONE=0, BLACK=30, BLUE=34, GREEN=32, CYAN=36, 
                  RED=31, PURPLE=35, YELLOW=33, WHITE=37};

  LoggerStream() 
  {
    colorMode_=true;
    color_=NONE;
    stream_=&std::cout;
    mute_=false;
  }
  void setStream(std::ostream* stream)
  {
    stream_=stream;
  }

  template <typename T>
  LoggerStream& operator << (const T& x)
  { 
    if (!mute_) 
      *stream_ << stringColor_ << x << stringReset_; 
    return *this;
  }

  LoggerStream& operator<<(std::ostream&(*f)(std::ostream&) )
  {
    if (!mute_)
      *stream_ << f;
    return *this;
  }

  void EnableColor()
  { colorMode_=true; update(); }

  void DisableColor()
  { colorMode_=false; update(); }

  void Mute()
  { mute_=true; }

  void UnMute()
  { mute_=false; }

  void SetColor(ColorType color)
  { color_=color; update(); }

  std::streamsize width() const
  { stream_->width(); }

  std::streamsize width (std::streamsize wide)
  { stream_->width(wide); }

private:

  void update()
  {
    if (!colorMode_  || color_==NONE)
    {
      stringColor_="";
      stringReset_="";
    }
    else
    {
      std::stringstream str;
      str << "\x1b[" << static_cast<unsigned int>(color_) << "m";
      stringColor_=str.str();
      stringReset_="\x1b[0m";
    }
  }


  std::ostream* stream_; 
  bool colorMode_;
  ColorType color_;
  bool mute_;
  std::string stringColor_;
  std::string stringReset_;

};




class Logger
{
public:

  static Logger* getInstance()
  {
    if (logger_==0) logger_ = new Logger;
    return logger_;
  }
  static void kill()
  {
    if (logger_!=0) delete logger_;
    logger_=0;
  }
  
  LoggerStream& debug()
  { return debug_; }
  LoggerStream& info()
  { return info_; }
  LoggerStream& warning(std::string fileName="",
                        signed int nLine=-1)
  {
    header(warning_,"WARNING",fileName,nLine);
    return warning_;
  }
  LoggerStream& error(std::string fileName="",
                      signed int nLine=-1)
  {
    header(error_,"ERROR",fileName,nLine);
    return error_;
  }
 
private:
  Logger() 
  {
    error_.SetColor(LoggerStream::RED);
    warning_.SetColor(LoggerStream::RED);
    info_.SetColor(LoggerStream::NONE);
    debug_.SetColor(LoggerStream::YELLOW);
  }

  ~Logger() {}
  static Logger* logger_;
  LoggerStream debug_;
  LoggerStream info_;
  LoggerStream warning_;
  LoggerStream error_;

  void header(LoggerStream& log,
              std::string logType,
              std::string fileName,
              signed int nLine)
  {
    log << logType;
    if (fileName!="" || nLine!=-1)
    {
      log << " ( ";
      log << fileName << " @ line=" << nLine;
      log << " )";
    }
    log << " : ";
  }

};


#endif
