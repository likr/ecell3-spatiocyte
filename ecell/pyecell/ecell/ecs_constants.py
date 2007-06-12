#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#       This file is part of the E-Cell System
#
#       Copyright (C) 1996-2007 Keio University
#       Copyright (C) 2005-2007 The Molecular Sciences Institute
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# E-Cell System is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# E-Cell System is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public
# License along with E-Cell System -- see the file COPYING.
# If not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
# 
#END_HEADER

#
# constants for ecell
#
from ecell._ecs import EntityType

# boolean constants
TRUE  = 1
FALSE = 0

# FullPN field numbers
TYPE       = 0
SYSTEMPATH = 1
ID         = 2
PROPERTY   = 3

# Entity type numbers
ENTITY     = EntityType.ENTITY
VARIABLE   = EntityType.VARIABLE
PROCESS    = EntityType.PROCESS
SYSTEM     = EntityType.SYSTEM

ENTITYTYPE_LIST = (
    EntityType.NONE,
    EntityType.ENTITY,
    EntityType.VARIABLE,
    EntityType.PROCESS,
    EntityType.SYSTEM,
    )

ENTITYTYPE_STRING_LIST = tuple([
    EntityType.getString( val ) for val in ENTITYTYPE_LIST ])

ENTITYTYPE_DICT = dict( zip( ENTITYTYPE_STRING_LIST, ENTITYTYPE_LIST ) )

# File extension

MODEL_FILE_EXTENSION = 'eml' 
SCRIPT_FILE_EXTENSION = 'ess' 


