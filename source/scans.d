/* Module for the scans class.
 *
 * Author: Jonathan Samson
 * Date: 04-08-2020
 */
module scans;

class mzXMLFile
/// Holds information relevant for an mzXML file, including a scan list.
/// Not fully implemented.
{
	real start_time;
	real end_time;
	string[] parent_file;
	string[] instrument_information;
	string[] data_processing;
	MSXScan[] scans;
}

/// Holds information relevant for a MS1 scan.
class Scan
{
	int scan_number;
	int centroided;
	uint level;
	string polarity;
	real retention_time;
	real[real] peaks;

	/* Debating whether these variables are necessary:
	string scanType;
	real total_ion_current;
	*/

	real get_peak_intensity(real my_peak)
	/* Gives the peak intensity of the set peak in the scan.
	 * Arguments:
	 *	my_peak - The peak of interest.
	 * Returns:
	 *	intensity - This scan's intensity of my_peak.
	 */
	{
		real* peak_intensity = (my_peak in peaks);
		real intensity = 0;
		if (peak_intensity !is null)
			intensity = *peak_intensity;
		return intensity;
	}

	void add_peak(real mz, real intensity)
	/* Either adds a peak or changes the intensity of a peak.
	 * Arguments:
	 *	mz - the mass to charge ratio of the new peak.
	 *	intensity - the intensity of the new peak.
	 */
	{
		peaks[mz] = intensity;
	}
}
unittest
{
	Scan test = new Scan;
	real[real] peaks = [
		51.46782684: 1460.6981201172,
		75.82749939: 1671.7169189453,
		75.86730194: 1605.3143310547,
		100.1144104: 1462.4990234375,
		101.5387802: 1490.517578125,
		107.7608643: 1808.1832275391,
		118.443428: 1619.8599853516,
		130.0875244: 37516.33203125,
		146.9610138: 1678.8117675781,
		171.1526642: 1760.8597412109,
		199.1815948: 35382.921875,
		243.1713562: 107272.828125,
		244.1736908: 8717.1875
	];
	test.level = 1;
	assert(test.level == 1);
	test.retention_time = 100.110;
	assert(test.retention_time == 100.110);
	test.peaks = peaks;
	assert(test.peaks == peaks);
	test.add_peak(56.12356, 5235.12359);
	peaks[56.12356] = 5235.12359;
	assert(test.get_peak_intensity(56.12356) == 5235.12359);
	assert(test.peaks == peaks);
	test.polarity = "-";
	assert(test.polarity == "-");
	assert(test.polarity != "+");
	test.scan_number = 1;
	assert(test.scan_number == 1);
	test.centroided = 1;
	assert(test.centroided == 1);
}

/// A Scan that includes parent data and collision energy.
class MSXScan : Scan
{
	Scan parent_scan;
	real parent_peak;
	float collision_energy;

	real get_parent_rt()
	/* Gives the retention time from the parent scan for this MSX.
	 * Returns:
	 *	this.parent_scan.get_rt() - The relevant parent rt.
	 */
	{
		return parent_scan.retention_time;
	}
}
unittest
{
	Scan parent = new Scan;
	parent.retention_time = 600.100;
	Scan notparent = new Scan;
	MSXScan test = new MSXScan;
	real[real] peaks = [
		51.46782684: 1460.6981201172,
		75.82749939: 1671.7169189453,
		75.86730194: 1605.3143310547,
		100.1144104: 1462.4990234375,
		101.5387802: 1490.517578125,
		107.7608643: 1808.1832275391,
		118.443428: 1619.8599853516,
		130.0875244: 37516.33203125,
		146.9610138: 1678.8117675781,
		171.1526642: 1760.8597412109,
		199.1815948: 35382.921875,
		243.1713562: 107272.828125,
		244.1736908: 8717.1875
	];
	test.level = 2;
	assert(test.level == 2);
	test.retention_time = 100.110;
	assert(test.retention_time == 100.110);
	test.peaks = peaks;
	assert(test.peaks == peaks);
	test.add_peak(56.12356, 5235.12359);
	peaks[56.12356] = 5235.12359;
	assert(test.get_peak_intensity(56.12356) == 5235.12359);
	assert(test.peaks == peaks);
	test.parent_peak = 244.1736908;
	assert(test.parent_peak == 244.1736908);
	test.parent_scan = parent;
	assert(test.parent_scan == parent);
	assert(test.parent_scan != notparent);
	assert(test.get_parent_rt() == 600.100);
	test.collision_energy = 35.0;
	assert(test.collision_energy == 35.0);
}
