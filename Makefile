
# --------------------------------------------------------------------------------
# 
# Copyrigtht 2012 Georgios Karagiannis
# 
# This file is part of AISRJ_Nmix.
# 
# AISRJ_Nmix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# AISRJ_Nmix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with AISRJ_Nmix.  If not, see <http://www.gnu.org/licenses/>.
# 
# --------------------------------------------------------------------------------

# References:
# 
# Karagiannis, G., & Andrieu, C. (2013). 
# Annealed importance sampling reversible jump MCMC algorithms. 
# Journal of Computational and Graphical Statistics, 22(3), 623-648.
# 
# Georgios Karagiannis
# School of Mathematics, University of Bristol
# University Walk, Bristol, BS8 1TW, UK
# Email : Georgios.Karagiannis@pnnl.gov
# Email (current): georgios-stats@gmail.com
# 
# Christophe Andrieu
# School of Mathematics, University of Bristol
# University Walk, Bristol, BS8 1TW, UK
# Email: C.Andrieu@bristol.ac.uk

FC = gfortran  #  choose a fortran compiler
FFLAGS = -O3   #  choose a flag

# Generate the executable

aisRJ.exe: algama.o rnguniform.o genprm.o rngnormal.o rnggamma.o rngbeta.o fixedmove_mod.o splitmergemove_mod.o aisrjnmix_pro.o
	$(FC) $(FFLAGS) -o aisrj.exe algama.o rnguniform.o genprm.o rngnormal.o rnggamma.o rngbeta.o fixedmove_mod.o splitmergemove_mod.o aisrjnmix_pro.o

algama.o: algama.f
	$(FC) $(FFLAGS) -c algama.f

rnguniform.o: rnguniform.f
	$(FC) $(FFLAGS) -c rnguniform.f
   
genprm.o: rnguniform.o genprm.f
	$(FC) $(FFLAGS) -c genprm.f
   
rngnormal.o: rnguniform.o rngnormal.f
	$(FC) $(FFLAGS) -c rngnormal.f
   
rnggamma.o: rnguniform.o rnggamma.f90
	$(FC) $(FFLAGS) -c rnggamma.f90

rngbeta.o: rnggamma.o rngbeta.f90
	$(FC) $(FFLAGS) -c rngbeta.f90
   
fixedmove_mod.mod: fixedmove_mod.o fixedmove_mod.f90 rngbeta.o rnggamma.o rngnormal.o genprm.o rnguniform.o algama.o
	$(FC) $(FFLAGS) -c fixedmove_mod.f90 
fixedmove_mod.o: fixedmove_mod.f90 rngbeta.o rnggamma.o rngnormal.o genprm.o rnguniform.o algama.o
	$(FC) $(FFLAGS) -c fixedmove_mod.f90  
   
splitmergemove_mod.mod: splitmergemove_mod.o splitmergemove_mod.f90 rngbeta.o rnggamma.o rngnormal.o genprm.o rnguniform.o algama.o
	$(FC) $(FFLAGS) -c splitmergemove_mod.f90 
splitmergemove_mod.o: splitmergemove_mod.f90 rngbeta.o rnggamma.o rngnormal.o genprm.o rnguniform.o algama.o
	$(FC) $(FFLAGS) -c splitmergemove_mod.f90
   
aisrjnmix_pro.o: fixedmove_mod.mod splitmergemove_mod.mod rngbeta.o rnggamma.o rngnormal.o genprm.o rnguniform.o algama.o aisrjnmix_pro.f90
	$(FC) $(FFLAGS) -c aisrjnmix_pro.f90

# Clean auxiliary files

clean:
	rm fixedmove_mod.mod fixedmove_mod.o splitmergemove_mod.mod splitmergemove_mod.o rngbeta.o rnggamma.o rngnormal.o genprm.o rnguniform.o algama.o aisrjnmix_pro.o


