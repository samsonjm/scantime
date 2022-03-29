/* Module for simulating TopN.
 *
 * Author: Jonathan Samson
 * Date: 28-03-2022
 */
module topn;
import mzxmlparser;
import scans;
import std.algorithm.searching;
import std.algorithm;

real[real] select_precursor_ions_topn(
	MSXScan[] scans,
	int n,
	real dew,
	real minimum_intensity,
	real mass_isolation_window,
	real full_scan_time)
{
/* Outputs list of [RT:mass] of precursors chosen for fragmentation
 * Arguments:
 *  scans - a MSXScan[] of only MS1 scans of the same polarity
 *  n - number of fragmentations per full scan
 *  dew - time between fragmentation selection of a single precursor ion
 *	minimum_intensity - minimum intensity to consider for fragmentaiton
 *  mass_isolation_window - m/z +/- to be considered the same peak
 *  full_scan_time - time between full scans when fragmenting
 * Returns:
 *  selected - [RT:mass] of precursor ions chosen
 * 
 */
	real[real] selected; // RT:M/Z  (FINAL list)
	real[real] peaks; // M/Z:intensity (all from current scan)
	real[real] top_peaks; // intensity:M/Z (current scan to add to final)
	real rt_offset = full_scan_time / (n + 1); // for fragmnetation time
	for(int scan_number = 0; scan_number<scans.length; ++scan_number)
	{
		real rt = scans[scan_number].retention_time;
		peaks = scans[scan_number].peaks;
		foreach(mz; peaks.byKey())
		{
			real new_intensity = peaks[mz];
			if(new_intensity < minimum_intensity)
			{
				continue;
			}
			if(top_peaks.length < n &&
			    scan_number == 0)
			{
				top_peaks[new_intensity] = mz;
				continue;
			}
			real[] intensities = top_peaks.keys;
			real lowest = 0;
			if(intensities.length > 0)
			{
				lowest = minElement(intensities);
			}
			if(new_intensity > lowest || top_peaks.length < n)
			{
				real[] selected_rts = selected.keys;
				int[] remove_list;
				for(int i=0; i<selected_rts.length; ++i)
				{
					if(selected_rts[i] < (rt - dew) ||
					   selected[selected_rts[i]] > (mz + mass_isolation_window) ||
					   selected[selected_rts[i]] < (mz - mass_isolation_window))
					{
						remove_list ~= i; 
					}
				}
				if(remove_list.length < selected_rts.length)
				{
					remove_list = [];
					continue;
				}
				if(top_peaks.length >= n)
				{
					top_peaks.remove(lowest);
				}
				top_peaks[new_intensity] = mz;
			}
		}
		if(top_peaks.length == 0)
		{
			continue;
		}
		auto sorted_intensities = (top_peaks.keys.sort.reverse);
		for(int i=0; i<top_peaks.length; ++i)
		{
			real adjusted_rt = rt + (rt_offset * (i + 1));
			selected[adjusted_rt] = top_peaks[sorted_intensities[i]];
		}
		top_peaks.clear();
	}
	return selected;
}
unittest
{
	MSXScan first = new MSXScan;
	first.level=1;
	first.retention_time = 1;
	real[real] first_peaks = [
		50.0: 1000000.0,
		55.0: 1000001.0,
		60.0: 1000002.0,
		65.0: 1000003.0,
		70.0: 1000004.0,
		75.0: 1000005.0,
		80.0: 1000006.0,
		85.0: 1000007.0,
		90.0: 1000008.0,
		95.0: 1000009.0,
		100.0: 1000010.0,
		105.0: 1000011.0,
		110.0: 1000012.0,
	];
	first.peaks = first_peaks;
	MSXScan second = new MSXScan;
	second.level=1;
	second.retention_time = 4;
	real[real] second_peaks = [
		50.0: 1000000.0,
		55.0: 1000001.0,
		60.0: 1000002.0,
		65.0: 1000003.0,
		70.0: 1000004.0,
		75.0: 1000005.0,
		80.0: 1000006.0,
		85.0: 1000007.0,
		90.0: 1000008.0,
		95.0: 1000009.0,
		100.0: 1000009.999,
		105.0: 1000011.001,
		110.0: 1000012.0,
	];
	second.peaks = second_peaks;
	MSXScan third = new MSXScan;
	third.level=1;
	third.retention_time = 5;
	real[real] third_peaks = [
		50.0: 1000000.0,
		55.0: 1000001.0,
		60.0: 1000002.0,
		65.0: 1000003.0,
		70.0: 1000004.0,
		75.0: 1000005.0,
		80.0: 1000006.0,
		85.0: 1000007.0,
		90.0: 1000008.0,
		95.0: 1000009.0,
		100.0: 1000010.0,
		105.0: 1000011.0,
		110.0: 1000012.0,
	];
	third.peaks = third_peaks;
	MSXScan fourth = new MSXScan;
	fourth.level=1;
	fourth.retention_time = 10;
	real[real] fourth_peaks = [
		50.0: 1000.0,
		55.0: 1001.0,
		60.0: 1002.0,
		65.0: 1003.0,
		70.0: 1004.0,
		75.0: 1005.0,
		80.0: 1006.0,
		85.0: 1007.0,
		90.0: 1008.0,
		95.0: 1009.0,
		100.0: 1010.0,
		105.0: 1011.0,
		110.0: 1012.0,
	];
	fourth.peaks = fourth_peaks;
	MSXScan[] all_scans = [first, second, third, fourth];
	int n = 4;
	int dew = 3;
	real min_intensity = 1000000.0;
	real mass_isolation_window = 0.01;
	int full_scan_time = 1;
	real[real] my_precursors = select_precursor_ions_topn(all_scans,
								  n,
								  dew,
								  min_intensity,
								  mass_isolation_window,
								  full_scan_time);
	real[] precursor_rts = my_precursors.keys;
	int num_per_scan = 0;
	foreach(rt; precursor_rts)
	{
		if(rt > 4 && rt < 5)
		{
			++num_per_scan;
		}
	}
	assert(num_per_scan == n);
	bool same_within_dew = false;
	foreach(rt, mz; my_precursors)
	{
		if(mz == 110.0 &&
		   rt > 2 &&
		   rt < 5)
		{
			same_within_dew = true;
		}

	}
	assert(!same_within_dew);
	foreach(rt; precursor_rts)
	{
		assert(rt < 10);
	}
	assert(my_precursors[1.2] == 110.0);
	assert(my_precursors[1.4] == 105.0);
	assert(my_precursors[1.6] == 100.0);
	assert(my_precursors[1.8] == 95.0);
	assert(my_precursors[4.2] < 1000008.5);
	assert(my_precursors[4.4] < 1000008.5);
	assert(my_precursors[4.6] < 1000008.5);
	assert(my_precursors[4.8] < 1000008.5);
}
