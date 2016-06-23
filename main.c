/*
    Tool to read data from an specific accelerometer (Phidget Prcision 003 High 
    Resoluation). The data can be stored into a wave file, but the main purpose
    is to apply a STFT and to create a CSV file per day.  
    
    It is written to run on Linux. It is only tested with Ubuntu and Debian.
	
    Copyright (C) 2015  Steffen Kühn / steffen.kuehn@em-sys-dev.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <glib.h>
#include <stdio.h>
#include <phidget21.h>
#include <unistd.h>
#include <sndfile.h>
#include <libgen.h>
#include <string.h>
#include <rfftw.h>
#include <math.h>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#define DEFAULT_OUTPUT_DIR "."
#define OUTPUT_MARKER "accel"
#define MAX_G 0.005 // 1.0 in wav is this value in g
#define DEFAULT_MAX_FREQ 150
#define DEFAULT_AVERAGE_INTERVAL_IN_SECONDS 10
#define PIPELINE_LEN 100

static SNDFILE* wavfile = 0;
static SF_INFO sfinfo = {0};
static char* output_dir = DEFAULT_OUTPUT_DIR;
static gboolean info_only = FALSE;
static int samplerate = 1000;
static int maxfreq = DEFAULT_MAX_FREQ;
static int avg_int_in_sec = DEFAULT_AVERAGE_INTERVAL_IN_SECONDS;
static int rbufi = 0;
static fftw_real*** inbuf = NULL;
static gboolean unproc[PIPELINE_LEN] = {0};
static int ibptr = 0;
static int aind = 0;
static fftw_real*** ampspec = NULL;
static gboolean max_instead_of_avg = FALSE;
static gboolean wav = FALSE;
static double moving_average[3] = {0};
static double tau = 10.0;// time in seconds in which the moving average is down to 0.5
static double avgconst = 0;

static GOptionEntry entries[] = {
	{
		"output-directory", 'd', 0, G_OPTION_ARG_FILENAME, &output_dir,
		"output dir, default: " DEFAULT_OUTPUT_DIR, NULL
	},
	{
		"info", 'i', 0, G_OPTION_ARG_NONE, &info_only,
		"show device info and terminate", NULL
	},
	{
		"average-interval", 'a', 0, G_OPTION_ARG_INT, &avg_int_in_sec,
		"averaging interval in seconds, default: " STR(DEFAULT_AVERAGE_INTERVAL_IN_SECONDS), NULL
	},
	{
		"max-frequency", 'm', 0, G_OPTION_ARG_INT, &maxfreq,
		"max. frequency in Hz, default: " STR(DEFAULT_MAX_FREQ), NULL
	},
	{
		"calcmax", 'M', 0, G_OPTION_ARG_NONE, &max_instead_of_avg,
		"calculate maximum instead of average", NULL
	},
	{
		"wav", 'w', 0, G_OPTION_ARG_NONE, &wav,
		"store wav file too", NULL
	},
	{ NULL}
};

static void calc_amplitude_spectrum(fftw_real* in, int N, fftw_real* amplitude_spectrum)
{
	fftw_real out[N];
	rfftw_plan p;
	int k;

	p = rfftw_create_plan(N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);

	rfftw_one(p, in, out);
	amplitude_spectrum[0] = sqrt(out[0] * out[0]);// DC component
	for (k = 1;k < (N + 1) / 2;++k) // (k < N/2 rounded up)
	{
		amplitude_spectrum[k] = sqrt(out[k] * out[k] + out[N - k] * out[N - k]);
	}
	if (N % 2 == 0) // N is even
	{
		amplitude_spectrum[N / 2] = sqrt(out[N / 2] * out[N / 2]);// Nyquist freq.
	}

	rfftw_destroy_plan(p);
}

static void free_spec_buffers(void)
{
	if (inbuf)
	{
		for (int j = 0;j < PIPELINE_LEN;j++)
		{
			for (int i = 0;i < 3;i++) g_free(inbuf[j][i]);
			g_free(inbuf[j]);
		}
		g_free(inbuf);
		inbuf = NULL;
	}

	if (ampspec)
	{
		for (int i = 0;i < 3;i++)
		{
			for (int j = 0;j < avg_int_in_sec;j++)
			{
				g_free(ampspec[i][j]);
			}
		}
		g_free(ampspec);
		ampspec = NULL;
	}
}

static void alloc_spec_buffers(void)
{
	free_spec_buffers();
	inbuf = g_new0(fftw_real**, PIPELINE_LEN);
	for (int j = 0;j < PIPELINE_LEN;j++)
	{
		inbuf[j] = g_new0(fftw_real*, 3);
		for (int i = 0;i < 3;i++)
		{
			inbuf[j][i] = g_new0(fftw_real, samplerate);
		}
	}
	ampspec = g_new0(fftw_real**, 3);
	for (int i = 0;i < 3;i++)
	{
		ampspec[i] = g_new0(fftw_real*, avg_int_in_sec);
		for (int j = 0;j < avg_int_in_sec;j++)
		{
			ampspec[i][j] = g_new0(fftw_real, samplerate / 2 + 1);
		}
	}
}

static void open_output(void)
{
	alloc_spec_buffers();
}

static gboolean does_file_exist(char* name)
{
	gboolean exist = FALSE;
	FILE* ofp = fopen(name, "r");
	if (ofp != NULL)
	{
		exist = TRUE;
		fclose(ofp);
	}
	return exist;
}

static gboolean csv_prepare(char* name)
{
	// check if the output file already exist
	char* open_mode = "w";
	if (does_file_exist(name)) open_mode = "a";

	// create output file
	FILE* ofp = fopen(name, open_mode);

	if (!ofp)
	{
		printf("ERROR: could not open/create output file: %s\n", name);
		return FALSE;
	}

	// create header
	if (strcmp(open_mode, "a"))
	{
		fprintf(ofp, "timestamp");
		for (int i = 0;i < maxfreq + 1;i++)
		{
			fprintf(ofp, ",%i Hz", i);
		}
		fprintf(ofp, "\n");
	}

	fclose(ofp);

	return TRUE;
}

static gboolean output_csv(int dim)
{
	time_t rawtime;
	struct tm * ti;
	char dims[2] = {0};

	time(&rawtime);
	ti = localtime(&rawtime);

	if (dim == 0)
	{
		dims[0] = 'x';
	}
	else if (dim == 1)
	{
		dims[0] = 'y';
	}
	else
	{
		dims[0] = 'z';
	}

	char* filename = NULL;
	filename = g_strdup_printf("%s/%4.4i-%2.2i-%2.2i_%s_%s.csv", output_dir,
		ti->tm_year + 1900, ti->tm_mon + 1, ti->tm_mday, dims, OUTPUT_MARKER);

	if (!csv_prepare(filename)) return FALSE;

	FILE* ofp = fopen(filename, "a");
	if (ofp == NULL)
	{
		printf("ERROR: could not reopen output file: %s\n", filename);
		return FALSE;
	}

	g_free(filename);

	fprintf(ofp, "%4.4i-%2.2i-%2.2i %2.2i:%2.2i:%2.2i",
		ti->tm_year + 1900, ti->tm_mon + 1, ti->tm_mday,
		ti->tm_hour, ti->tm_min, ti->tm_sec);

	for (int k = 0;k < maxfreq + 1;k++)
	{
		float v = 0;
		if (max_instead_of_avg)
		{
			for (int j = 0;j < avg_int_in_sec;j++)
			{
				if ((ampspec[dim][j][k] > v) || (j == 0))
				{
					v = ampspec[dim][j][k];
				}
				v /= (samplerate / 1000.0);
			}
		}
		else
		{
			for (int j = 0;j < avg_int_in_sec;j++)
			{
				v += ampspec[dim][j][k];
			}
			v /= (avg_int_in_sec * samplerate / 1000.0);// unit is mg
		}

		fprintf(ofp, ",%f", v);
	}

	fprintf(ofp, "\n");
	fclose(ofp);

	return TRUE;
}

static void process(void)
{
	for (int k = 0;k < PIPELINE_LEN;k++)
	{
		int ptr = ibptr + k + (PIPELINE_LEN / 10);
		while ((ptr >= PIPELINE_LEN) && (ptr >= 0)) ptr = ptr - PIPELINE_LEN;

		if (unproc[ptr])
		{
			if (ptr >= PIPELINE_LEN) ptr = 0;

			for (int i = 0;i < 3;i++)
			{
				calc_amplitude_spectrum(inbuf[ptr][i], samplerate, ampspec[i][aind]);
			}
			aind++;

			if (aind == avg_int_in_sec)
			{
				aind = 0;

				for (int i = 0;i < 3;i++)
				{
					output_csv(i);
				}
			}

			unproc[ptr] = FALSE;
		}
	}
}

static void close_wav(void)
{
	if (wavfile != 0)
	{
		sf_close(wavfile);
		wavfile = 0;
	}
}

static void write_wav(double data[3], gboolean check)
{
	if (check)
	{
		time_t rawtime;
		time(&rawtime);
		struct tm* ti = localtime(&rawtime);

		// create file name
		char* filename = NULL;
		filename = g_strdup_printf("%s/%4.4i-%2.2i-%2.2i_%s.wav", output_dir,
			ti->tm_year + 1900, ti->tm_mon + 1, ti->tm_mday, OUTPUT_MARKER);

		if (!does_file_exist(filename))
		{
			// close old file
			close_wav();

			// file does not exist. create it and open it		
			sfinfo.channels = 3;
			sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_DOUBLE;
			sfinfo.samplerate = samplerate;
			wavfile = sf_open(filename, SFM_WRITE, &sfinfo);

			// initialise moving average
			for (int i = 0;i < 3;i++) moving_average[i] = data[i];
		}

		if (wavfile == 0)
		{
			// file exists but it is not open
			wavfile = sf_open(filename, SFM_RDWR, &sfinfo);
			sf_seek(wavfile, 0, SEEK_END);

			// initialise moving average
			for (int i = 0;i < 3;i++) moving_average[i] = data[i];
		}
	}

	// calculate moving average
	for (int i = 0;i < 3;i++)
	{
		moving_average[i] = avgconst * moving_average[i] + (1.0 - avgconst) * data[i];
		data[i] -= moving_average[i];
	}

	if (wavfile != 0)
	{
		sf_write_double(wavfile, data, 3);
	}
}

// callback that will run when data arrive
int CCONV SpatialDataHandler(CPhidgetSpatialHandle spatial, void *userptr, 
	CPhidgetSpatial_SpatialEventDataHandle *data, int count)
{
	for (int k = 0;k < count;k++)
	{
		if (wav)
		{
			double buf[3] = {
				data[k]->acceleration[0] / MAX_G, // x
				data[k]->acceleration[1] / MAX_G, // y
				data[k]->acceleration[2] / MAX_G  // z 	
			};		
			write_wav(buf, k == 0);
		}

		for (int i = 0;i < 3;i++)
		{
			double v = data[k]->acceleration[i];
			inbuf[ibptr][i][rbufi] = v;
		}
		rbufi++;

		if (rbufi == samplerate)
		{
			rbufi = 0;
			unproc[ibptr] = TRUE;
			ibptr++;
			if (ibptr >= PIPELINE_LEN) ibptr = 0;

			if (unproc[ibptr])
			{
				printf("Realtime error!\n");
			}
		}
	}
	return 0;
}

static void close_output(void)
{
	close_wav();
}

// callback that will run if the sensor is attached to the computer
int CCONV AttachHandler(CPhidgetHandle spatial, void *userptr)
{
	int serialNo;
	CPhidget_getSerialNumber(spatial, &serialNo);
	printf("Spatial %10d attached!\n", serialNo);

	return 0;
}

// callback that will run if the sensor is detached from the computer
int CCONV DetachHandler(CPhidgetHandle spatial, void *userptr)
{
	int serialNo;
	CPhidget_getSerialNumber(spatial, &serialNo);
	printf("Spatial %10d detached!\n", serialNo);

	return 0;
}

// callback that will run if the sensor generates an error
int CCONV ErrorHandler(CPhidgetHandle spatial, void *userptr, int ErrorCode, const char *unknown)
{
	printf("Error handled. %d - %s\n", ErrorCode, unknown);
	return 0;
}

static int display_properties(CPhidgetHandle phid)
{
	int serialNo, version;
	const char* ptr;
	int numAccelAxes, numGyroAxes, numCompassAxes, dataRateMax, dataRateMin;

	CPhidget_getDeviceType(phid, &ptr);
	CPhidget_getSerialNumber(phid, &serialNo);
	CPhidget_getDeviceVersion(phid, &version);
	CPhidgetSpatial_getAccelerationAxisCount((CPhidgetSpatialHandle)phid, &numAccelAxes);
	CPhidgetSpatial_getGyroAxisCount((CPhidgetSpatialHandle)phid, &numGyroAxes);
	CPhidgetSpatial_getCompassAxisCount((CPhidgetSpatialHandle)phid, &numCompassAxes);
	CPhidgetSpatial_getDataRateMax((CPhidgetSpatialHandle)phid, &dataRateMax);
	CPhidgetSpatial_getDataRateMin((CPhidgetSpatialHandle)phid, &dataRateMin);

	printf("%s\n", ptr);
	printf("Serial Number: %10d\nVersion: %8d\n", serialNo, version);
	printf("Number of Accel Axes: %i\n", numAccelAxes);
	printf("Number of Gyro Axes: %i\n", numGyroAxes);
	printf("Number of Compass Axes: %i\n", numCompassAxes);
	printf("datarate> Max: %d  Min: %d\n", dataRateMax, dataRateMin);

	return 0;
}

static gboolean controlloop(void)
{
	int result;
	const char *err;
	CPhidgetSpatialHandle spatial = 0;

	avgconst = pow(2.0, -1.0 / (tau * samplerate));

	// create the spatial object
	CPhidgetSpatial_create(&spatial);

	// set the event handlers
	CPhidget_set_OnAttach_Handler((CPhidgetHandle)spatial, AttachHandler, NULL);
	CPhidget_set_OnDetach_Handler((CPhidgetHandle)spatial, DetachHandler, NULL);
	CPhidget_set_OnError_Handler((CPhidgetHandle)spatial, ErrorHandler, NULL);

	// open the spatial object for device connections
	CPhidget_open((CPhidgetHandle)spatial, -1);

	// get the program to wait for a spatial device to be attached
	printf("Waiting for spatial to be attached.... \n");
	for (;;)
	{
		if ((result = CPhidget_waitForAttachment((CPhidgetHandle)spatial, 10000)))
		{
			CPhidget_getErrorDescription(result, &err);
			printf("Problem waiting for attachment: %s\n", err);
			// wait five seconds for the next run
			usleep(5000000);
		}
		else
		{
			break;
		}
	}

	if (info_only)
	{
		// show properties of the attached sensor
		display_properties((CPhidgetHandle)spatial);
	}
	else
	{
		open_output();

		// register data callback
		CPhidgetSpatial_set_OnSpatialData_Handler(spatial, SpatialDataHandler, NULL);

		// set data rate for the spatial events
		CPhidgetSpatial_setDataRate(spatial, (int)(1000 / samplerate));

		for (;;)
		{
			process();
			usleep(2000);
		}
	}

	close_output();

	printf("Closing...\n");
	CPhidget_close((CPhidgetHandle)spatial);
	CPhidget_delete((CPhidgetHandle)spatial);

	return TRUE;
}

int main(int argc, char *argv[])
{
	gboolean res = TRUE;
	GOptionContext *context = NULL;

	context = g_option_context_new("");
	g_option_context_set_summary(context, "reads acceleration data from a \"Phidget Spatial 003 High Resolution\"-sensor");
	g_option_context_add_main_entries(context, entries, NULL);

	if (!g_option_context_parse(context, &argc, &argv, NULL))
	{
		res = FALSE;
		printf("ERROR: invalid options\n");
	}
	else if (!controlloop())
	{
		res = FALSE;
		printf("ERROR: problem in mainloop\n");
	}

	// return code 0 means everything was OK
	return (!res);
}
