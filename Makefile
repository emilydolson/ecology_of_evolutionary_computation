# Project-specific settings
PROJECT := evo_comp_ecology
EMP_DIR := ../Empirical/source

# Flags to use regardless of compiler
CFLAGS_all := -Wall -Wno-unused-function -std=c++1z -I$(EMP_DIR)/ -I../cppitertools

# Native compiler information
CXX := clang++-6.0
CXX_nat := $(CXX)
CFLAGS_nat := -O3 -DNDEBUG $(CFLAGS_all)
CFLAGS_nat_debug := -g $(CFLAGS_all)

# Emscripten compiler information
CXX_web := emcc
OFLAGS_web_all := -s TOTAL_MEMORY=67108864 --js-library $(EMP_DIR)/web/library_emp.js --js-library $(EMP_DIR)/web/d3/library_d3.js -s EXPORTED_FUNCTIONS="['_main', '_empCppCallback']" -s DISABLE_EXCEPTION_CATCHING=1 -s NO_EXIT_RUNTIME=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['cwrap', 'stringToUTF8', 'writeStringToMemory']"#--embed-file configs
OFLAGS_web := -Oz -DNDEBUG
OFLAGS_web_debug := -Oz -pedantic -Wno-dollar-in-identifier-extension

CFLAGS_web := $(CFLAGS_all) $(OFLAGS_web) $(OFLAGS_web_all)
CFLAGS_web_debug := $(CFLAGS_all) $(OFLAGS_web_debug) $(OFLAGS_web_all)

default: $(PROJECT)
native: $(PROJECT)
web: $(PROJECT).js
all: $(PROJECT) $(PROJECT).js

debug:	CFLAGS_nat := $(CFLAGS_nat_debug)
debug:	$(PROJECT)

debug-web:	CFLAGS_web := $(CFLAGS_web_debug)
debug-web:	$(PROJECT).js

web-debug:	debug-web

$(PROJECT):	source/native/$(PROJECT).cc
	$(CXX_nat) $(CFLAGS_nat) source/native/$(PROJECT).cc -o $(PROJECT)
	@echo To build the web version use: make web

$(PROJECT).js: source/web/$(PROJECT)-web.cc
	$(CXX_web) $(CFLAGS_web) source/web/$(PROJECT)-web.cc -o web/$(PROJECT).js

interaction_networks:	interaction_networks.cc
	$(CXX_nat) $(CFLAGS_nat) interaction_networks.cc -o interaction_networks

debug-interaction_networks:	interaction_networks.cc
	$(CXX_nat) $(CFLAGS_nat_debug) interaction_networks.cc -o interaction_networks

interaction_networks.js: interaction_networks-web.cc
	$(CXX_web) $(CFLAGS_web) interaction_networks-web.cc -o web/interaction_networks.js

debug-interaction_networks.js: interaction_networks-web.cc
	$(CXX_web) $(CFLAGS_web_debug) interaction_networks-web.cc -o web/interaction_networks.js

test: test_interaction_networks.cc
	$(CXX_nat) $(CFLAGS_nat_debug) test_interaction_networks.cc -o test
	./test


clean:
	rm -f $(PROJECT) web/$(PROJECT).js *.js.map *~ source/*.o

# Debugging information
print-%: ; @echo '$(subst ','\'',$*=$($*))'
