import std.stdio;
import scans;
import mscosine;
import calculatescantimes;
import mzxmlparser;
import std.getopt;
import std.regex;
import std.typecons;

void main(string[] args)
{
	string input_file;
	string summary_file_name;
    auto helpInformation = getopt(
                args,
                "input|i", "The input file in .mgl or .mzxml format, must be MS/MS of N>=2",
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
    string file_contents = read_file(input_file);
    auto file_extension = ctRegex!(`\.(\w*)$`);
    MSXScan[] my_scans;
    switch (input_file.matchFirst(file_extension)[1])
    {
        default:
        {
            throw new Exception("Invalid input file extension.");
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
	real[] scan_times = calculate_scan_times(my_scans, summary_file_name);
	real full_scan_time = scan_times[0];
	real fragmentation_scan_time = scan_times[1];
	writeln("Times:");
	writeln("Full Scan: ", full_scan_time);
	writeln("Fragmentation Scan: ", fragmentation_scan_time);
}
