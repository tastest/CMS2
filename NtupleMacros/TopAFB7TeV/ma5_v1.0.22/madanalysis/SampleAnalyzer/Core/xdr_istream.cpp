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


#include "Core/xdr_istream.h"

// -----------------------------------------------------------------------------
// overloading >> operator for std::string
// -----------------------------------------------------------------------------
xdr_istream& xdr_istream::operator>>(std::string &s)
{
  if (eof()) return (*this);
	UInt_t len;
	(*this)>>len;

  if (eof()) return (*this);

	char line[len];
	sb_->sgetn(line, len);
	s=std::string(line,line+len);
	
	size_t pad = 4-len&3;

	char dummy[pad];
	sb_->sgetn(dummy,pad);
	return *this;
}

// -----------------------------------------------------------------------------
// overloading >> operator for UInt_t
// -----------------------------------------------------------------------------
xdr_istream& xdr_istream::operator>>(UInt_t &v)
{
  if (eof()) return (*this);

	v=0;
	for(int i=0;i<32; i+=8)
  {
		v<<=8;
		v += static_cast<UInt_t>(sb_->sbumpc());
	}
	return *this;
}

// -----------------------------------------------------------------------------
// overloading >> operator for Int_t
// -----------------------------------------------------------------------------
xdr_istream& xdr_istream::operator>>(Int_t &v)
{
  if (eof()) return (*this);

	UInt_t _v;
	(*this)>>_v;
	v=static_cast<Int_t>(_v);
	return (*this);
}


// -----------------------------------------------------------------------------
// overloading >> operator for ULong_t
// -----------------------------------------------------------------------------
xdr_istream& xdr_istream::operator>>(ULong_t &v)
{
  if (eof()) return (*this);

	v=0;
	for(int i=0;i<64; i+=8){
		v<<=8;
		v += static_cast<ULong_t>(sb_->sbumpc());
	}
	return *this;
}


// -----------------------------------------------------------------------------
// overloading >> operator for Long_t
// -----------------------------------------------------------------------------
xdr_istream& xdr_istream::operator>>(Long_t &v)
{
  if (eof()) return (*this);

	ULong_t _v;
	(*this)>>_v;
	v=static_cast<Long_t>(_v);
	return (*this);
}


// -----------------------------------------------------------------------------
// overloading >> operator for Float_t
// -----------------------------------------------------------------------------
xdr_istream& xdr_istream::operator>>(Float_t &v)
{
  if (eof()) return (*this);

	UInt_t n;
	(*this)>>n;
	Float_t* vp = reinterpret_cast<Float_t*>(&n);
	v=*vp;
	return *this;
}


// -----------------------------------------------------------------------------
// overloading >> operator for Double_t
// -----------------------------------------------------------------------------
xdr_istream& xdr_istream::operator>>(Double_t &v)
{
  if (eof()) return (*this);

	ULong_t n;
	(*this)>>n;
	Double_t* vp = reinterpret_cast<Double_t*>(&n);
	v=*vp;
	return *this;
}
