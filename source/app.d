import std.stdio;
import scans;
import mscosine;
import calculatescantimes;
import mzxmlparser;
import std.getopt;
import std.typecons;

void main(string[] args)
{
	string input_file;
	string summary_file_name;
    auto helpInformation = getopt(
                args,
                "input|i", "The input files in .mgl or .mzxml format, must be MS/MS of N>=2, separated by a single comma (,)",
                &input_file,
				"summary_file|s", "The file to save the statistical summary in.  This file will be overwritten", 
				&summary_file_name);
    if(helpInformation.helpWanted)
    {
        defaultGetoptFormatter(
                    stdout.lockingTextWriter(),
                    "Outputs time for full and fragmentation scans",
                    helpInformation.options,
                    "  %*s\t%*s%*s%s\n");
        return;
    }
	real[] scan_times = calculate_scan_times(input_file, summary_file_name);
	real full_scan_time = scan_times[0];
	real fragmentation_scan_time = scan_times[1];
	writeln("Times:");
	writeln("Full Scan: ", full_scan_time);
	writeln("Fragmentation Scan: ", fragmentation_scan_time);
}
