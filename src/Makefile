# Requires:
#   BOOST libraries: boost.org
#   Gnu scientific libraries (GSL): https://www.gnu.org/software/gsl
#   HDF5 library: http://www.hdfgroup.org/HDF5/
#   Semi-Markov library: https://github.com/afidd/Semi-Markov
#   mcrand library: https://github.com/afidd/mcrand
#
# mcrand is a pain. It exists so that runs with different seeds will
# be reliably random, but it depends on RNGSSELIB by Barash and Shchur,
# which you have to retrieve from Elsevier to build it. If you feel lucky with
# the Boost random generator, then go to sir_exp.hpp and change
# RandGen to use that one, and you can take mcrand out of this Makefile.

BOOST=/usr/local/boost_1_54_0mt
# Different Boost installations have different suffixes.
# If there is no suffix, use "BOOSTVARIANT=".
BOOSTVARIANT=
SEMIMARKOV=/usr/local/include/semimarkov-0.1
HDF5=/usr/local/hdf5-1.8.11
MCRAND=$(HOME)/Documents/mcrand
#FIGTREE=$(HOME)/Documents/figtree
#FIGTREE_LIBS=-L$(FIGTREE)/lib -lann_figtree_version -lfigtree

-include local.mk
# -DSMVHIDELOG -pg
OPT=-std=c++11 -DHAVE_CONFIG_H -fPIC -g -O2 -DSMVHIDELOG
INCLUDES=-I$(SEMIMARKOV) -I. -I$(BOOST)/include -I$(HDF5)/include -I$(MCRAND)/include
LIBS=-L$(BOOST)/lib -L$(HDF5)/lib -L$(MCRAND)/lib \
    -lboost_unit_test_framework$(BOOSTVARIANT) \
	-lboost_log_setup$(BOOSTVARIANT) -lboost_log$(BOOSTVARIANT) \
	-lboost_chrono$(BOOSTVARIANT) -lboost_thread$(BOOSTVARIANT) \
	-lboost_date_time$(BOOSTVARIANT) -lboost_filesystem$(BOOSTVARIANT) \
	-lboost_program_options$(BOOSTVARIANT) -lboost_random$(BOOSTVARIANT) \
	-lboost_system$(BOOSTVARIANT) -lboost_timer$(BOOSTVARIANT) \
	-lhdf5 -lhdf5_hl -lgsl -lgslcblas \
	-lmcrand -lpthread

all: seir rider individual ensemble_sum
.PHONY: all

seir: seir_exp.o main.o hdf_file.o trajectory.o pen.o rider_enums.o \
		place_transition.o model_options.o
	$(CXX) $(OPT) -o seir seir_exp.o main.o  hdf_file.o trajectory.o pen.o \
	rider_enums.o place_transition.o model_options.o $(LIBS)

seir_exp.o: seir_exp.cpp seir_exp.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

main.o: main.cpp feedlot_version.hpp seir_exp.hpp hdf_file.hpp
	$(CXX) main.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o main.o

rider: rider.o ridermain.o hdf_file.o trajectory.o pen.o rider_enums.o \
		place_transition.o model_options.o
	$(CXX) $(OPT) -fPIC -o rider rider.o ridermain.o  hdf_file.o \
	trajectory.o pen.o rider_enums.o place_transition.o model_options.o $(LIBS)

rider.o: rider.cpp rider.hpp
	$(CXX) rider.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o rider.o

ridermain.o: ridermain.cpp feedlot_version.hpp rider.hpp hdf_file.hpp
	$(CXX) ridermain.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o ridermain.o


individual: individual.cpp fmdv.hpp feedlot_places.hpp feedlot_individual_transitions.hpp
	$(CXX) individual.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-o individual $(LIBS)

ensemble_sum: ensemble_sum.o hdf_file.o
	$(CXX) ensemble_sum.o hdf_file.o $(INCLUDES) $(OPT) \
	-o $@ $(LIBS)

ensemble_sum.o: ensemble_sum.cpp hdf_file.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

hdf_file.o: hdf_file.cpp hdf_file.hpp trajectory.hpp parameter.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

trajectory.o: trajectory.cpp trajectory.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

pen.o: pen.cpp pen.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

rider_enums.o: rider_enums.cpp rider_enums.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

place_transition.o: place_transition.cpp place_transition.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

model_options.o: model_options.cpp model_options.hpp
	$(CXX) $< $(OPT) $(INCLUDES) -c -o $@

feedlot_version.hpp: Makefile
	python getgit.py

clean:
	rm -f *.o sirexp together rider individual ensemble_sum sirdemo_version.hpp
