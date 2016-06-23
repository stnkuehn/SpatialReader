#!/usr/bin/python

# Plot CSV files which was generated with "spatialreader".
#
# Copyright (C) 2015  Steffen Kuehn / steffen.kuehn@em-sys-dev.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import csv
import argparse
import Image

import matplotlib.pyplot as plt
import matplotlib.cm as cm


def append_Hz(v):
	return "%i Hz" % v


def get_hour(str):
	return int(str[11:13])


def shorten_timestamp(str):
	year = str[0:4]
	month = str[5:7]
	day = str[8:10]
	hour = str[11:13]
	minute = str[14:16]
	return "%s-%s-%s %s:%s" % (year, month, day, hour, minute)


# returns 1 for time2 > time1, 0 for time2 == time1 and -1 for time2 < time1
def compare_timestamps(time2, time1):
	hour2 = int(time2[0:2])
	minute2 = int(time2[3:5])
	hour1 = int(time1[0:2])
	minute1 = int(time1[3:5])
	if hour2 > hour1:
		return 1
	elif hour2 < hour1:
		return -1
	if minute2 > minute1:
		return 1
	elif minute2 < minute1:
		return -1
	return 0


def read_csvfile(name, mintime, maxtime, minfreq, maxfreq, dboffset, scale):
	f = open(name, 'r')
	reader = csv.DictReader(f)
	frequencies = None
	date = None
	values = []
	timestamps_index = []
	timestamps_labels = []
	ohour = None
	index = 0
	timecut = True

	if compare_timestamps(mintime, maxtime) == 0:
		timecut = False

	lines = iter(reader)
	linenbr = 0

	while True:
		linenbr = linenbr + 1

		try:
			row = next(lines)
		except Exception as e:
			if len(str(e)) == 0:
				# end of file
				print('file %s ends at line %d' % (name, linenbr))
				break
			else:
				print('error in file %s, line %d: %s' % (name, linenbr, str(e)))
				pass
		else:
			# extract head line
			if frequencies is None:
				frequencies = row.keys()
				frequencies.remove('timestamp')
				frequencies = map(int, map(lambda c: c[:-3], frequencies))
				frequencies.sort()
				frequencies = filter(lambda a: a >= int(minfreq), frequencies)
				frequencies = filter(lambda a: a <= int(maxfreq), frequencies)
				freqs = map(append_Hz, frequencies)

			timestamp = row['timestamp']
			tswy = timestamp[11:16]
			cmin = compare_timestamps(mintime, tswy)  # should not be 1
			cmax = compare_timestamps(maxtime, tswy)  # should not be -1
			if timecut and (cmin == 1 or cmax == -1):
				continue

			hour = get_hour(timestamp)
			if ohour is None or hour is not ohour:
				ohour = hour
				timestamps_index.append(index)
				timestamps_labels.append(shorten_timestamp(timestamp))

			# matrix with sorted values
			line = []
			try:
				for key in freqs:
					v = float(row[key])
					v -= float(dboffset)
					v *= float(scale)
					line.append(v)
				values.append(line)
			except Exception as e:
				print('error in file %s, line %d: %s' % (name, linenbr, str(e)))
				pass

			index = index + 1

	f.close()
	values = zip(*values)
	return (frequencies, timestamps_index, timestamps_labels, values)


def plot_csvfile(filename, tdpi, mintime, maxtime, odir, minfreq,
                 maxfreq, aspectratio, freqdist, downscale, dboffset, maxdb, legendtext,
                 scale):
	print("process file %s" % filename)
	(filewoext, ext) = os.path.splitext(filename)
	(path, fname) = os.path.split(filewoext)
	if odir is not None:
		path = odir
	if path == '':
		path = '.'
	imagename = path + '/' + fname + '.jpg'
	(frequencies, timestamps_index, timestamps_labels, values) = read_csvfile(filename, mintime, maxtime, minfreq,
	                                                                          maxfreq, dboffset, scale)
	if len(values) > 0:
		pixelx = len(values[0])
		pixely = float(aspectratio) * len(values)
		fig = plt.figure(figsize=(pixelx / 100, pixely / 100))
		img = plt.imshow(values, cmap=cm.jet, clim=(0, float(maxdb)), aspect=float(aspectratio), origin='lower')
		ax = plt.gca()
		ax.set_xticks(timestamps_index)
		ax.set_xticklabels(timestamps_labels)
		filter_freqs = filter(lambda a: a % int(freqdist) == 0, frequencies)
		filter_inds = []
		for filter_freq in filter_freqs:
			filter_inds.append(frequencies.index(filter_freq))
		ax.set_yticks(filter_inds)
		ax.set_yticklabels(filter_freqs)
		cbar = plt.colorbar(img, pad=0.01)
		cbar.ax.get_yaxis().labelpad = 15
		cbar.ax.set_ylabel(legendtext, rotation=270)
		plt.ylabel('frequency [Hz]', fontdict={'fontsize': 15})
		plt.xlabel('time of day', fontdict={'fontsize': 15})
		plt.savefig(imagename, dpi=float(tdpi), bbox_inches='tight')
		plt.close(fig)
		# downsampling
		factor = float(downscale)
		if factor < 1.0 and factor > 0.0:
			img_org = Image.open(imagename)
			width_org, height_org = img_org.size
			img_ds = img_org.resize((int(factor * width_org), int(factor * height_org)), Image.ANTIALIAS)
			img_ds.save(imagename)


def get_csv_files_in_directory(dir):
	res = []
	for path, subdirs, files in os.walk(dir):
		for name in files:
			(root, ext) = os.path.splitext(name)
			if ext.lower() == '.csv':
				res.append(os.path.join(path, name))
	res.sort()
	return res


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--file', help="CSV file with dB values")
	parser.add_argument('-d', '--dir', help="process all csv files in this directory")
	parser.add_argument('-o', '--odir', help="output directory")
	parser.add_argument('-r', '--dpi', help="output image dpi value", default='200')
	parser.add_argument('-a', '--mintime', help="start time, default: 00:00", default='00:00')
	parser.add_argument('-b', '--maxtime', help="stop time, default: 24:00", default='24:00')
	parser.add_argument('-A', '--minfreq', help="minimal frequency", default='0')
	parser.add_argument('-B', '--maxfreq', help="maximal frequency", default='1000')
	parser.add_argument('-R', '--aspectratio', help="aspect ratio", default='2.5')
	parser.add_argument('-D', '--freqdist', help="frequency distance in scala", default='10')
	parser.add_argument('-s', '--downscale', help="downscale factor 0..1", default='1')
	parser.add_argument('-O', '--dboffset', help="dB threshold", default='0')
	parser.add_argument('-m', '--maxdb', help="upper limit for dB", default='40')
	parser.add_argument('-l', '--legendtext', help="text for legend", default='')
	parser.add_argument('-S', '--scale', help="scalefactor (for unit change)", default='1.0')
	args = parser.parse_args()

	if args.dir is not None:
		files = get_csv_files_in_directory(args.dir)
	else:
		files = [args.file]

	for file in files:
		plot_csvfile(file, args.dpi, args.mintime, args.maxtime, args.odir, args.minfreq, args.maxfreq,
		             args.aspectratio, args.freqdist, args.downscale, args.dboffset, args.maxdb,
		             args.legendtext, args.scale)


if __name__ == "__main__":
	main()
