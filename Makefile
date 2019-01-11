EXE	= ePoPE
CC	= gcc
LL	= gcc
VERSION	= 2.0
WARN	=  -Wall -Werror
#CCFLAGS	= $(WARN) -O3 -O1
CCFLAGS	= $(WARN) -g -std="c99"
LIBS	= -lm
SRCS	= main.c read.c tree.c calc.c newicktolist2testbaumok.c
OBJS	= $(foreach file, $(SRCS:.c=.o),$(file))

all: $(EXE)
$(EXE): $(OBJS)
	$(LL) $(CCFLAGS) $(OBJS) $(LIBS) -o $(EXE)
$(OBJS): %.o: %.c
	 $(CC) -c $(CCFLAGS) $< -o $@

clean:
	rm -f $(OBJS) core $(EXE) *~*
