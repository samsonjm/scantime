import std.stdio;
import scans;
import mscosine;
import topn;
import mzxmlparser;
import std.getopt;
import std.regex;

void main(string[] args)
{
	string input_file;
	int n;
	real dew;
	int min_intensity;
	real mass_iso;
	real full_scan_time;
	bool filter_c13_isotopologs;
    auto helpInformation = getopt(
                args,
                "input|i", "The input file in .mgl or .mzxml format",
                &input_file,
                "N value|n", "The N of TopN", &n,
                "DEW|d", "The DEW, in seconds", &dew,
				"min_intensity|m", "The minimum intensity for fragmentation",
				&min_intensity,
				"mass_iso_window|w", "The mass isolation width", &mass_iso,
				"full_scan_time|s", "The time in seconds for a full scan", &full_scan_time);
		 		//"filter_c13_isotopologs|f", "'true' to filter C13 isotopologs", &filter_c13_isotopologs);
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
        case "mgl":
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
	real[real] selected_precursors = select_precursor_ions_topn(my_scans,
			n,
			dew,
			min_intensity,
			mass_iso,
			full_scan_time);
			//filter_c13_isotopologs);
	writeln("RT\tM/Z");
	foreach(rt, mz; selected_precursors)
	{
		writefln("%s\t%s", rt, mz);
	}
}
