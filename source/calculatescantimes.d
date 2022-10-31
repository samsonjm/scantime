/* Module for simulating TopN.
 *
 * Author: Jonathan Samson
 * Date: 12-0r73-2022
 */
module calculatescantimes;
import mzxmlparser;
import scans;
import mscosine;
import std.algorithm.iteration : uniq;
import std.algorithm.mutation : copy;
import dstats;
import std.stdio;
import std.regex;
import std.array : array, split;
import std.algorithm : map;
import std.range : zip, chain, repeat;

import ggplotd.aes : aes;
import ggplotd.axes : xaxisLabel, yaxisLabel, xaxisRange, yaxisRange;
import ggplotd.geom : geomBox;
import ggplotd.ggplotd : GGPlotD, putIn, Margins, title;

real[] calculate_scan_times(
		string input_files,
		string summary_file_name
		)
{
/* Outputs list of [0] full and [1] fragmentation scan times averaged over files
 *  Arguments:
 *   input_files - comma-separated input files
 *   summary_file_name - location to save summary file
 * Returns:
 *  scan_times - a list of [0] full and [1]  fragmentation scan times
 */
	string[] separate_files = input_files.split(",");
	real[] full_scan_times;
	real[] fragmentation_scan_times;
	auto file_extension = ctRegex!(`\.(\w*)$`);
	MSXScan[] my_scans;
	foreach(file; separate_files)
	{
		string file_contents = read_file(file);
		switch (file.matchFirst(file_extension)[1])
		{
			default:
			{
				throw new Exception("Invalid input file  extension.");
			}
			case "mgf":
			{
				my_scans = mgf_parser(file_contents);
				break;
			}
			case "mzXML":
			{
				my_scans = parse_mzxml(file_contents);
				break;
			}
		}
		real[] scan_times = calculate_separate_scan_times(my_scans, summary_file_name);
		full_scan_times ~= scan_times[0];
		fragmentation_scan_times ~= scan_times[1];
	}
	real full_scan_time = mean(full_scan_times);
	real fragmentation_scan_time = mean(fragmentation_scan_times);
	return [full_scan_time, fragmentation_scan_time];
}
unittest
{
}

real[] calculate_separate_scan_times(
		MSXScan[] scans,
		string summary_file_name)
{
/* Outputs list of [0] full and [1] fragmentation scan times
 * Arguments:
 *  scans - a MSXScan[] of only MS1 scans of the same polarity
 * Returns:
 *  scan_times - a list of [0] full and [1] fragmentation scan times
 */
	File summary_file = File(summary_file_name, "w");
	real rt;
	int level;
	real[] full_scan_times;
	real[] frag_scan_times;
	real[] double_frag_scan_times;
	real[] triple_frag_scan_times;
	real[] quadruple_frag_scan_times;
	real[] last_frag_scan_times;
	for(int scan_number = 0; scan_number < scans.length; ++scan_number)
	{
		rt = scans[scan_number].retention_time;	
		level = scans[scan_number].level;
		if(scan_number < scans.length - 1)
		{
			real rt_diff = scans[scan_number+1].retention_time - rt;
			if(level == 1)
			{
				full_scan_times ~= rt_diff;
			}
			else if(scans[scan_number+1].level == 2)
			{
				frag_scan_times ~= rt_diff;
			}
			else if(scans[scan_number+1].level == 1)
			{
				last_frag_scan_times ~= rt_diff;
			}
		}
	}
	summary_file.writeln("Summary of full scan times:");
	summary_file.writeln(summary(full_scan_times));
	summary_file.writeln("\nLinear regression of full scan times:");
	summary_file.writeln(linearRegress(full_scan_times, repeat(mean(full_scan_times))));
	summary_file.writeln("\nSummary of fragmentation scan times:");
	summary_file.writeln(summary(frag_scan_times));
	summary_file.writeln("\nLinear regression of fragmentation scan times:");
	summary_file.writeln(linearRegress(frag_scan_times, repeat(1)));
	summary_file.writeln("\nSummary of last-fragmentation-before-a-full-scan scan times:");
	summary_file.writeln(summary(last_frag_scan_times));
	summary_file.writeln("\nLinear regression of last fragmentation scan times:");
	summary_file.writeln(linearRegress(last_frag_scan_times, repeat(1)));
	real[] full_scan_uniq = full_scan_times;
	real[] frag_scan_uniq = frag_scan_times;
	real[] last_frag_scan_uniq = last_frag_scan_times;
	full_scan_uniq.length -= full_scan_uniq.uniq().copy(full_scan_uniq).length;
	frag_scan_uniq.length -= frag_scan_uniq.uniq().copy(frag_scan_uniq).length;
	last_frag_scan_uniq.length -= last_frag_scan_uniq.uniq().copy(last_frag_scan_uniq).length;
	auto plot_array = full_scan_uniq ~ frag_scan_uniq ~  last_frag_scan_uniq;
	auto cols = "full scan".repeat(full_scan_uniq.length)
		        .chain("frag scan".repeat(frag_scan_uniq.length))
				.chain("last frag scan".repeat(last_frag_scan_uniq.length));
	auto gg = plot_array.zip(cols)
							 .map!((a) => aes!("x","colour","fill","label")(a[0], a[1], 0.45, a[1]))
							 .geomBox
							 .putIn(GGPlotD());
	gg.put(yaxisLabel("Time (s)"));
	gg.put(title("Scan Times - Unique values"));
	gg.save("/home/samsonjm/Projects/ScanTime/boxplot.svg");
	real frag_scan = sum(frag_scan_times) / frag_scan_times.length;
	real full_scan = (sum(full_scan_times) / full_scan_times.length);
	real last_frag = sum(last_frag_scan_times) / last_frag_scan_times.length;
	full_scan += last_frag;
	full_scan -= frag_scan;
	return [full_scan, frag_scan];
}
unittest
{
	MSXScan first = new MSXScan;
	first.level=1;
	first.retention_time = 1;
	MSXScan second = new MSXScan;
	second.level=2;
	second.retention_time=1.75;
	MSXScan third = new MSXScan;
	third.level=2;
	third.retention_time=2.25;
	MSXScan fourth = new MSXScan;
	fourth.level=2;
	fourth.retention_time=2.75;
	MSXScan fifth = new MSXScan;
	fifth.level=2;
	fifth.retention_time=3.25;
	MSXScan sixth = new MSXScan;
	sixth.level=2;
	sixth.retention_time=3.75;
	MSXScan seventh = new MSXScan;
	seventh.level=1;
	seventh.retention_time=4.5;
	MSXScan eighth = new MSXScan;
	eighth.level=2;
	eighth.retention_time=5.25;
	MSXScan ninth = new MSXScan;
	ninth.level=2;
	ninth.retention_time=5.75;
	MSXScan tenth = new MSXScan;
	tenth.level=2;
	tenth.retention_time=6.25;
	MSXScan eleventh = new MSXScan;
	eleventh.level=2;
	eleventh.retention_time=6.75;
	MSXScan twelvth = new MSXScan;
	twelvth.level=2;
	twelvth.retention_time=7.25;
	MSXScan thirteenth = new MSXScan;
	thirteenth.level=1;
	thirteenth.retention_time=8;
	MSXScan fourteenth = new MSXScan;
	fourteenth.level=9;
	fourteenth.retention_time=8.75;
	MSXScan fifteenth = new MSXScan;
	fifteenth.level=2;
	fifteenth.retention_time=9.25;
	MSXScan sixteenth = new MSXScan;
	sixteenth.level=2;
	sixteenth.retention_time=9.75;
	MSXScan seventeenth = new MSXScan;
	seventeenth.level=2;
	seventeenth.retention_time=10.25;
	MSXScan eighteenth = new MSXScan;
	eighteenth.level=2;
	eighteenth.retention_time=10.75;
	MSXScan[] all_scans = [first, second, third, fourth, fifth, sixth,
						   seventh, eighth, ninth, tenth, eleventh, twelvth,
						   thirteenth, fourteenth, fifteenth, sixteenth, seventeenth, eighteenth];
	assert(calculate_separate_scan_times(all_scans, "/home/samsonjm/Projects/ScanTime/summary.txt") == [1, 0.5]);
}
