# To build vcf_parser:
#
#   make
#
# run with
#
#   ./build/vcf_parser

D_COMPILER=ldc2

LDMD=ldmd2

DUB_LIBS =

DFLAGS = -wi -I./source $(DUB_INCLUDE)
RPATH  =
LIBS   =
SRC    = $(wildcard source/vcf_parser/*.d  source/test/*.d)
IR     = $(wildcard source/vcf_parser/*.ll source/test/*.ll)
BC     = $(wildcard source/vcf_parser/*.bc source/test/*.bc)
OBJ    = $(SRC:.d=.o)
OUT    = build/vcf_parser

debug: DFLAGS += -O0 -g -d-debug $(RPATH) -link-debuglib $(BACKEND_FLAG) -unittest
release: DFLAGS += -O -release $(RPATH)

.PHONY:test

all: debug

build-setup:
	mkdir -p build/

ifeq ($(FORCE_DUPLICATE),1)
  DFLAGS += -d-version=FORCE_DUPLICATE
endif


default debug release profile getIR getBC gperf: $(OUT)

# ---- Compile step
%.o: %.d
	$(D_COMPILER) -lib $(DFLAGS) -c $< -od=$(dir $@) $(BACKEND_FLAG)

# ---- Link step
$(OUT): build-setup $(OBJ)
	$(D_COMPILER) -of=build/vcf_parser $(DFLAGS)  $(OBJ) $(LIBS   =) $(DUB_LIBS) $(BACKEND_FLAG)

test:
	chmod 755 build/vcf_parser
	./run_tests.sh

debug-strip: debug

clean:
	rm -rf build/*
	rm -f $(OBJ) $(OUT) trace.{def,log}
