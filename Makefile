TARGET = spatialreader

OBJECTS = main.o

PKGS = glib-2.0

CFLAGS = -std=gnu99 -Wall -funsigned-char `pkg-config --cflags $(PKGS)` -DSTATIC=static
LDFLAGS = `pkg-config --libs $(PKGS)` -lm -lrfftw -lfftw -lphidget21 -lsndfile

ifdef DEBUG
	CFLAGS += -ggdb -O0
else
	CFLAGS += -O2 -Werror
endif

all: $(TARGET)

# link
$(TARGET): $(OBJECTS)
	 $(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)
clean:
	-rm -f $(OBJECTS) $(TARGET)

