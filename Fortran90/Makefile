COMP = gfortran
FLAGS = -Wall -O3
LIBS = -L/opt/local/lib/ -lslatec -llapack -lblas -latlas

########
# MAIN #
########

HF:
	$(COMP) $(FLAGS) $(LIBS) input.f90 output.f90 constants.f90 factorial.f90 gaussian.f90 la.f90 overlap.f90 multipole.f90 kinetic.f90 nuclear.f90 core.f90 electronic.f90 density.f90 energy.f90 init.f90 DIIS.f90 RHF.f90 UHF.f90 force.f90  HF.f90 -o HF.x


################################

dynamics: BOMD

molecules: H2Test HeHTest H2OTest

tests: factorialTest gaussianTest

BOMD:
	$(COMP) $(FLAGS) $(LIBS) input.f90 output.f90 constants.f90 factorial.f90 gaussian.f90 la.f90 overlap.f90 multipole.f90 kinetic.f90 nuclear.f90 core.f90 electronic.f90 density.f90 energy.f90 init.f90 DIIS.f90 RHF.f90 UHF.f90 force.f90 dynamics.f90  tests/BOMD_test.f90 -o BOMD.x

factorialTest:
	$(COMP) $(FLAGS) $(LIBS) factorial.f90 tests/factorial_test.f90 -o factorialTest.x

gaussianTest:
	$(COMP) $(FLAGS) $(LIBS) factorial.f90 gaussian.f90 tests/gaussian_test.f90 -o gaussianTest.x

H2Test:
	$(COMP) $(FLAGS) $(LIBS) constants.f90 factorial.f90 gaussian.f90 overlap.f90 kinetic.f90 nuclear.f90 core.f90 electronic.f90 utils.f90 density.f90 energy.f90 RHF.f90 input.f90 tests/H2_test.f90 -o H2Test.x

HeHTest:
	$(COMP) $(FLAGS) $(LIBS) constants.f90 factorial.f90 gaussian.f90 overlap.f90 kinetic.f90 nuclear.f90 core.f90 electronic.f90 utils.f90 density.f90 energy.f90 RHF.f90 input.f90 tests/HeH_test.f90 -o HeHTest.x

H2OTest:
	$(COMP) $(FLAGS) $(LIBS) constants.f90 factorial.f90 gaussian.f90 overlap.f90 kinetic.f90 nuclear.f90 core.f90 electronic.f90 utils.f90 density.f90 energy.f90 RHF.f90 input.f90 tests/H2O_test.f90 -o H2OTest.x

clean:
	rm -f *.mod

cleanall:
	rm -f *.mod *.x
