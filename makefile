# Set compiler's options depending on the architecture
ifeq ($(ARCHI),i86)
  CC     = clang
#  CFLAGS = -g -c
  CFLAGS = -O3 -c -Wuninitialized -Wunused -Winline -Wshadow \
           -fexpensive-optimizations -funroll-loops
  LDFLAGS= -g 
endif
ifeq ($(ARCHI),ppc)
#  CC      = /opt/ibmcmp/vac/6.0/bin/xlc
#  CFLAGS  = -c -O3 -qarch=auto
#  LDFLAGS = -s
  CC     = gcc
  CFLAGS = -O3 -s -c \
	   -Wuninitialized -Wunused -Winline -Wshadow \
           -fexpensive-optimizations -funroll-loops
#CFLAGS = -g -c
#  LDFLAGS=
endif
ifeq ($(ARCHI),win32)
  CC    = gcc
  CFLAGS= -c -O3 -mno-cygwin -Wuninitialized -Wunused -Winline -Wshadow
  LDFLAGS=-v -s -mno-cygwin
endif

# working dirs
EXEDIR = build
SRCDIR = sources
OBJDIR = objects/$(ARCHI)
ARCDIR = archives
DIRDIR = objects $(EXEDIR) $(OBJDIR) $(ARCDIR)
INCDIR = 
LDLDIR = 
VPATH  = $(SRCDIR)

# objects list
src    = $(wildcard $(SRCDIR)/*.c)
header = $(wildcard $(SRCDIR)/*.h)
objs   = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(src:.c=.o))
prog   = mshdist

#.SILENT:

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(OPT64) $(INCDIR) $(CFLAGS) $< -o $@

$(EXEDIR)/$(prog):$(DIRDIR) $(objs)
	echo "#define COMPIL " '"' `date` '"' > $(SRCDIR)/compil.date
	$(CC) -c $(CFLAGS) $(INCDIR) $(SRCDIR)/mshdist.c -o $(OBJDIR)/mshdist.o
	$(CC) $(LDFLAGS) $(OPT64) $(LDLDIR) $(objs) -o $@ -lm

$(objs):$(header)

$(DIRDIR):
	@[ -d $@ ] || mkdir $@

clean:
	-rm $(objs) $(EXEDIR)/$(prog)

tar:$(DIRDIR)
	tar czf $(ARCDIR)/$(prog).`date +"%Y.%m.%d"`.tgz sources makefile

target: $(EXEDIR)/$(prog)
