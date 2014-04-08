gnomish
=======
a PERL scripting language extension for simulation of plant genetic stystems.


Gnomish is an engine for simulating genetic crossing at the trait level in plants.  

The application is mostly complete and functional.
However, there are some loose ends in the documentation and some incomplete features. 

The tutorial and examples should give the reader a fairly good idea of how to
write simulations with Gnomish.  Generally speaking, methods and functions 
in the various Gnomish packages take named parameters, which are fairly self explanatory.  
When in doubt, it is best to consult the code itself.

This package is free software and distributed under the GNU General Public License
as published by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

If you make any changes or updates to this package, the authors would
appreciate a copy of your changes, so that they can be incorporated
into future releases of Gnomish.

Authors:

Hai Pham, hpham42@gmail.com
Diane Mather, diane.mather@adelaide.edu.au
Nick Tinker, nick.tinker@agr.gc.ca  

November 27, 2005


-------------------------------------------------------------------------------------

NOTES:

This is a preliminary release of Gnomish.  The tutorials in this package are 
currently only designed to be executed in the 'gnomish' subdirectory.  
Please do not move them anywhere else.

For example, if you choose to unzip the gnomish.zip file in your C: drive, 
then the gnomish subdirectory should be at 'C:\gnomish'.

To execute the tutorial simulation 'hardy-weinberg':

- Open a command line shell (click Start->Run, enter 'cmd' at the prompt)
- type 'cd \gnomish' in the command line shell
- type 'perl hardy-weinberg.pl'

Any new scripts can be created and executed from the \gnomish subdirectory 
as well.


INSTALLATION:

There is a Perl package for Windows 2000 and Windows XP in the 
\gnomish\perl subdirectory.  Please refer to the README.txt file there for 
installation instructions.


Documentation:

There are some documentation for the various Gnomish modules and the 
tutorials in \gnomish\docs.  These are preliminary documentation, and will 
need refinement.  Please make comments about them, particularly where they 
can be improved upon.


Database interface:

The Gnomish database interface will require the installation of the Perl 
modules DBI and DBD from the standard Perl module library.  I will need to 
study this in somewhat more detail to create installation instructions for 
these modules.  There is an empty Access database file in \gnomish\database 
for examination.


To do:

- The Gnomish::Evalulation module is still very preliminary and will need 
to be visited as we spend more time with simulations and decide on the the 
commonly used evaluation methods and selection methods.  If we find 
ourselves custom coding the same evaluation/selection method repeatedly, it 
should probably go into the standard Gnomish evaluation library.

- Database access and installation still needs some work.

- Tutorial 2 leaves some open questions (these are noted within the text of 
the tutorial).  We will need to revisit these once we understand how the 
simulation engine will be used.

- User interface is a big open question.  In order to make Gnomish useful 
to the non-programming audience, we will definately need a user interface.  
Interfaces are also useful if we just want a quick answer to a problem or 
the ability to quickly setup and run a simple simulation.

