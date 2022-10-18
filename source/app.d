import std.stdio;
import scans;
import mscosine;
import topn;
import mzxmlparser;
import std.getopt;
import std.regex;
import std.typecons;

void main(string[] args)
{
	string input_file;
	int n;
	real dew;
	int min_intensity;
	real mass_iso;
	real total_scan_time;
	bool filter_c13_isotopologues;
	int max_c13_in_isotopologues;
	int max_charge;
	real lag_time;
    auto helpInformation = getopt(
                args,
                "input", "The input file in .mgf or .mzxml format - full scan only",
                &input_file,
                "N value|n", "The N of TopN", &n,
                "DEW|d", "The DEW, in seconds", &dew,
				"min_intensity|m", "The minimum intensity for fragmentation",
				&min_intensity,
				"total_scan_time|s", "The time between full scans in seconds", &total_scan_time,
		 		"filter_c13_isotopologues|f", "'true' to filter C13 isotopologues", &filter_c13_isotopologues,
				"mass_iso_window|w", "The mass isolation width in ppm - only used by isotopologue filter", &mass_iso,
				"max_c13_in_isotopologues|i", "Maximum number of C13 isotopologues in a peak to filter (default=4)", &max_c13_in_isotopologues,
				"max_charge|c", "The maximum expected charge of the ions", &max_charge,
				"lag_time|l", "The time between detecting a M/Z and its first fragmentation", &lag_time);
    if(helpInformation.helpWanted)
    {
        defaultGetoptFormatter(
                    stdout.lockingTextWriter(),
                    "Runs in silico TopN precursor ion selection",
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
	Tuple!(real, real)[real] selected_precursors = select_precursor_ions_topn(my_scans,
			n,
			dew,
			min_intensity,
			mass_iso,
			total_scan_time,
			filter_c13_isotopologues,
			max_c13_in_isotopologues,
			max_charge,
			lag_time);
	writeln("RT\tM/Z\tIntensity");
	foreach(rt, tup; selected_precursors)
	{
		real mz = tup[0];
		real intensity = tup[1];
		writefln("%s\t%s\t%s", rt, mz, intensity);
	}
}
