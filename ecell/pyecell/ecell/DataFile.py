#!/usr/bin/env python2

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#        This file is part of E-CELL Session Monitor package
#
#                Copyright (C) 1996-2002 Keio University
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# E-CELL is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# E-CELL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public
# License along with E-CELL -- see the file COPYING.
# If not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
# 
#END_HEADER
#
# Design: Kouichi Takahashi <shafi@e-cell.org>
# Design and Programming: Masahiro Sugimoto <sugi@bioinformatics.org> at
# E-CELL Project, Lab. for Bioinformatics, Keio University.
#

from string import *

# ------------------------------------------------------------------
# DataFile (This is abstract class)
#   - manages one file object
#   - has file name
# ------------------------------------------------------------------
class DataFile:

	# ------------------------------------------------------------------
	# Constructor
	#
	# return  -> None
	# ------------------------------------------------------------------
	def __init__(self):

		self.theFileName = ''

	# end of __init__


	# ------------------------------------------------------------------
	# setData (abstract method)
	#
	# aData(object) : data object
	#
	# If this method doesn't be implemented in sub class, 
	# then throws NotImplementedError
	# ------------------------------------------------------------------
	#def setData(self, aData):

	#	import inspect
	#	caller = inspect.getouterframes(inspect.currentframe())[0][3]
	#	raise NotImplementedError(caller + 'must be implemented in subclass')

	# end of setData


	# ------------------------------------------------------------------
	# setFileName
	#
	# aFileName(string)  : a file name 
	#
	# return -> None
	# ------------------------------------------------------------------
	def setFileName(self, aFileName):

		if( type(aFileName) != type('') ):
			raise TypeError("Error : aFileName is not list matrix.")

		if( len(aFileName) == 0 ):
			print 'Warning: %s, the length of filename is 0' %__name__

		self.theFileName = aFileName

	# end of setFileName


	# ------------------------------------------------------------------
	# getFileName
	#
	# return -> the file name (string)
	# ------------------------------------------------------------------
	def getFileName(self):

		return self.theFileName

	# end of theFileName


	# ------------------------------------------------------------------
	# theData
	#
	# return -> the data(string list)
	# ------------------------------------------------------------------
	#def theData(self):

	#	import inspect
	#	caller = inspect.getouterframes(inspect.currentframe())[0][3]
	#	raise NotImplementedError(caller + 'must be implemented in subclass')

	# end of theData


	# ------------------------------------------------------------------
	# save ( abstract )
	#
	# If this method doesn't be implemented in sub class, 
	# then throws NotImplementedError
	# ------------------------------------------------------------------
	def save(self):

		import inspect
		caller = inspect.getouterframes(inspect.currentframe())[0][3]
		raise NotImplementedError(caller + 'must be implemented in subclass')

	# end of save


	# ------------------------------------------------------------------
	# saveWithFileName
	#
	# aFileName(string)  : a file name 
	#
	# return -> None
	# This method is throwable exception.
	# ------------------------------------------------------------------
	def saveWithFileName(self, aFileName):
	
		self.setFileName(aFileName)
		self.save()

	# end of saveWithFileName


	# ------------------------------------------------------------------
	# load
	#
	# If this method doesn't be implemented in sub class, 
	# then throws NotImplementedError
	# ------------------------------------------------------------------
	def load(self):

		import inspect
		caller = inspect.getouterframes(inspect.currentframe())[0][3]
		raise NotImplementedError(caller + 'must be implemented in subclass')

	# end of load


	# ------------------------------------------------------------------
	# loadWithFileName
	#
	# aFileName(string)  : a file name 
	#
	# return -> None
	# This method is throwable exception.
	# ------------------------------------------------------------------
	def loadWithFileName(self, aFileName):

		self.setFileName(aFileName)
		self.load()

	# end of loadWithFileName


# end of DataFile

if __name__ == "__main__":

	class SubClass1(DataFile):
		def setData(self, aData):
			print "setData"
		def theData(self):
			print "theData"
		def save(self):
			print "save"
		def load(self):
			print "load"

	class SubClass2(DataFile):
		def setData(self, aData):
			print "setData"
		def theData(self):
			print "theData"
		def save(self):
			print "save"
		def load(self):
			print "load"

	def main():
		sub = SubClass1()
		sub.setData('hoge')
		sub.theData()
		sub.save()
		sub.load()

		sub = SubClass2()
		sub.setData('hoge')
		sub.theData()
		sub.save()
		sub.load()

		file = open('hoge','w')
		file.close()

	main()

