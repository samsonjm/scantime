/* Tools to parse mzXML files into an MSXScan[].
 * 
 * Author: Jonathan Samson
 * Date: 04-08-2020
 */
module mzxmlparser;
import scans;
import std.bitmanip;
import std.conv;
import std.base64;
import std.stdio;
import std.math;
import std.exception;
import std.algorithm;
import dxml.parser;

real[real] decode_mzxml_string(
		string encoded, 
		string compression="none", 
		int precision=32)
{
/* Decodes the mzXML string that represents the scan peaks.
 * Arguments:
 *	encoded - The encoded string.
 *	compression - The type of compression (only zlib accepted).
 *	precision - The precision used to encode the string.
 * Returns:
 *	peak_list - the m/z:intensity list of all peaks in the sacn.
 */
	ubyte[] decoded = Base64.decode(encoded);
	real mz;
	real intensity;
	real[real] peak_list;
	int byte_size = precision / 8;
	enforce(compression == "none" || compression == "zlib",
			"Invalid compression type.");
	enforce(precision == 64 || precision == 32,
			"Invalid precision.");
	if (compression=="zlib")
	{
		import std.zlib;
		decoded = cast(ubyte[]) uncompress(decoded);
	}
	for(int i = 1; i<=decoded.length/byte_size; ++i)
	{
		float readable;
		if (precision == 64)
		{
			ubyte[8] next_value = decoded[byte_size*(i-1)..
						     byte_size*i];
			readable = bigEndianToNative!double(next_value);
		}
		else // precision = 32
		{
			ubyte[4] next_value = decoded[byte_size*(i-1)..
						     byte_size*i];
			readable = bigEndianToNative!float(next_value);
		}
		if(i % 2 == 0)
		{
			intensity = readable.to!real;
			peak_list[mz] = intensity;
		}
		else
			mz = readable.to!real;
	}
	return peak_list;
}
unittest
{
	import std.algorithm;
	import std.math;
	real[real] answer = [
		51.4678:	1460.7,
		75.8275:	1671.72,
		75.8673:	1605.31,
		100.114:	1462.5,
		101.539:	1490.52,
		107.761:	1808.18,
		118.443:	1619.86,
		130.088:	37516.3,
		146.961:	1678.81,
		171.153:	1760.86,
		199.182:	35382.9,
		243.171:	107273,
		244.174:	8717.19
	];
	string line = "Qk3fDkS2lldCl6euRND28UKXvA9EyKoPQsg6lES2z/hCyxPbRLp" ~
		"QkELXhZBE4gXdQuzjCUTKe4VDAhZoRxKMVUMS9gVE0dn6QysnFUTcG4ND" ~
		"Ry59Rwo27ENzK95H0YRqQ3Qsd0YINMA=";
	real[real] function_test = decode_mzxml_string(line);
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));
	line = "eJwBaACX/0JN3w5EtpZXQpenrkTQ9vFCl7wPRMiqD0LIOpREts/4QssT20" ~
		"S6UJBC14WQROIF3ULs4wlEynuFQwIWaEcSjFVDEvYFRNHZ+kMrJxVE3Bu" ~
		"DQ0cufUcKNuxDcyveR9GEakN0LHdGCDTAJ+wubA==";
	function_test = decode_mzxml_string(line, "zlib");
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));
	line = "QEm74cAAAABAltLK4AAAAEBS9PXAAAAAQJoe3iAAAABAUveB4AAAAECZFU" ~
		"HgAAAAQFkHUoAAAABAltn/AAAAAEBZYntgAAAAQJdKEgAAAABAWvCyAAA" ~
		"AAECcQLugAAAAQF2cYSAAAABAmU9woAAAAEBgQs0AAAAAQOJRiqAAAABA" ~
		"Yl7AoAAAAECaOz9AAAAAQGVk4qAAAABAm4NwYAAAAEBo5c+gAAAAQOFG3" ~
		"YAAAABAbmV7wAAAAED6MI1AAAAAQG6FjuAAAABAwQaYAAAAAA==";
	function_test = decode_mzxml_string(line, "none", 64);
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));
	line = "eJxz8Nz98AADA4PDtEunHoDooC9fwfxZcvcUwPzvjWDxmaKOYDqSPagBrP" ~
		"7mfwYwP6k6AURP9xIC86M+bALTcxx2LwDRsXMSwebM9C8A8xOczoLlHwV" ~
		"2gflJcQfA9CxrewcQnZryCMyf3VwANjfj6Xkw/6HbXbC9eanVYPf9MugF" ~
		"q89r7QO76yDbDJC5AD9eO64=";
	function_test = decode_mzxml_string(line, "zlib", 64);
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));
	assertThrown(decode_mzxml_string(line, "7z", 64));
	assertThrown(decode_mzxml_string(line, "none", 5));
}

string read_file(string name_of_file)
/* Reads the file into a string.
 * Arguments:
 *      file_stream - The name of the file to read.
 *
 * Returns:
 *      file_contents - The contents of the file.
 */
{
    string file_contents = "";
    try
    {
        auto file = File(name_of_file, "r");
        string line;
        while ((line = file.readln()) !is null)
        {
            file_contents ~= line;
        }
        file.close();
    }
    catch(ErrnoException e)
    {
        writeln("Invalid file name");
    }
    return file_contents;
}

MSXScan[] parse_mzxml(string contents)
/* Parses the contents of an .mzXML file into a list of Scan objects.
 * Arguments:
 *	contents - the contents of a .mzXML file.
 * Returns:
 *	scans - a list of Scan objects populated by contents.
 *
 * This parser uses dxml to parse the string.
 */
{
	auto range = parseXML(contents);
	range.skipToPath("scan");
	MSXScan[] scans;
	MSXScan current_scan;
	while(range.empty == false)
	{
		if (range.front.type == EntityType.elementStart)
		{
			auto attr = range.front.attributes;
			switch (range.front.name)
			{
				case "scan":
				{
					current_scan = new MSXScan;
					string retention_time;
					attr.getAttrs(
							"num", 
							&current_scan.scan_number,
							"centroided", 
							&current_scan.centroided,
							"msLevel", 
							&current_scan.level,
							"polarity", 
							&current_scan.polarity,
							"retentionTime", 
							&retention_time,
							"collisionEnergy", 
							&current_scan.collision_energy
							);
					current_scan.retention_time = 
						retention_time[2..$-1].to!real;
					break;
				}
				case "precursorMz":
				{
					int parent_scan;
					attr.getAttrs(
							"precursorScanNum", 
							&parent_scan);
					current_scan.parent_scan = 
						scans[parent_scan-1];
					break;
				}
				case "peaks":
				{
					string compression_type;
					int precision;
					//string byte_order;
					attr.getAttrs(
							"compressionType", 
							&compression_type,
							"precision", 
							&precision,
							//"byteOrder", 
							//&byte_order;
							);
					range.popFront();
					string encoded_read = range.front.text;
					current_scan.peaks = 
						decode_mzxml_string(
							encoded_read,
							compression_type,
							precision);
					break;
				}
				default:
					break;
			}
		}
		else if (range.front.type == EntityType.elementEnd)
		{
			string name  = range.front.name;
			if (name == "scan")
				scans ~= current_scan;	
			else if (name == "msRun")
				break;
		}
		range.popFront();
	}
	return scans;
}
unittest
{
	string scans = read_file("testfiles/example.mzXML");
	MSXScan[] parsed = parse_mzxml(scans);
	assert(parsed[0].scan_number == 1);
	assert(parsed[0].centroided == 1);
	assert(parsed[0].level == 1);
	assert(parsed[0].polarity == "-");
	assert(approxEqual(parsed[0].retention_time, 0.457129));
	assert(approxEqual(parsed[2].retention_time, 647.132));
	assert(parsed[0].collision_energy.isNaN);
	assert(parsed[2].collision_energy == 35.0);
	assert(parsed[0].parent_scan is null);
	assert(parsed[2].parent_scan == parsed[1]);
	real[real] peaks = [
 		51.4678:    1460.7,
        75.8275:    1671.72,
        75.8673:    1605.31,
        100.114:    1462.5,
        101.539:    1490.52,
        107.761:    1808.18,
        118.443:    1619.86,
        130.088:    37516.3,
        146.961:    1678.81,
        171.153:    1760.86,
        199.182:    35382.9,
        243.171:    107273,
        244.174:    8717.19
	];
	assert(approxEqual(parsed[2].peaks.keys.sort, peaks.keys.sort));
	assert(approxEqual(parsed[2].peaks.values.sort, 
				peaks.values.sort)); assert(parsed[2].level == 2);
	assert(approxEqual(parsed[2].get_peak_intensity(
					parsed[2].peaks.keys.sort[8]
					), 
				1678.81));
	assert(parsed[2].scan_number == 3);
}
