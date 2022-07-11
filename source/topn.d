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

bool is_c13_isotopologue(
		real current_mz,
		real former_mz,
		real mass_isolation_window,
		int charge=4,
		int max_isotopologues=1)
{
/* Returns true if the selected mz is a C13 isotopologue
 * Arguments:
 *	current_mz - The M/Z to check
 *	former_mz - a single M/Z to check against
 *	mass_isolation_window - M/Z +/- to be considered the same peak
 *	max_charge - The maximum expected charge
 *  max_isotopologues - The maximum expected 13C in a peak
 * Returns:
 *	isotopologue - boolean, TRUE if a C13 isotopologue
 * Only returns TRUE if the M/Z is a C13 isotopologue.  Returns FALSE
 * if the M/Z is the monoisotopic peak or if there are no other peaks
 * within the C13 isotopologue range.
 */
	real min_limit;
	real max_limit;
	real theoretical_monoisotopologue_mz;
	min_limit = former_mz - mass_isolation_window;
	max_limit = former_mz + mass_isolation_window;
	for(int num_c13 = 1; num_c13 < max_isotopologues + 1; ++num_c13)
	{
		for(int i = 1; i < charge + 1; ++i)
		{
			theoretical_monoisotopologue_mz = current_mz - 
				(num_c13 * (1.00335484 / i));
			if(theoretical_monoisotopologue_mz > min_limit &&
			   theoretical_monoisotopologue_mz < max_limit)
			{
				return true;
			}
		}
	}
	return false;
}
unittest
{
	real my_peak_mzs = 160;
	real mass_iso_window = 0.00000001;
	assert(!is_c13_isotopologue(160, // Doesn't detect self FALSE
								my_peak_mzs, 
								mass_iso_window, 
								1));
	assert(!is_c13_isotopologue(150, // No close-by peaks, FALSE
							   my_peak_mzs,
							   mass_iso_window,
							   1));
	assert(!is_c13_isotopologue(160, // Monoisotopic peak, FALSE
							   my_peak_mzs,
							   mass_iso_window,
							   1));
// M/Z outside of isotopologue range should be FALSE
	assert(!is_c13_isotopologue(161.00335486, // Just above C13 window, FALSE
							   my_peak_mzs,
							   mass_iso_window,
							   1));
	assert(!is_c13_isotopologue(161.00335482, // Just below C13 window, FALSE
							   my_peak_mzs,
							   mass_iso_window,
							   1));
// M/Z within isotopologue range should be TRUE
	assert(is_c13_isotopologue(161.003354831, // Just inside C13 window low end, TRUE
							   my_peak_mzs,
							   mass_iso_window,
							   1));
	assert(is_c13_isotopologue(161.00335485, // Just inside C13 window high end, TRUE
							   my_peak_mzs,
							   mass_iso_window,
							   1));
// Charge of 2 should work for 1/2 the window
	assert(is_c13_isotopologue(160.50167742, // Inside window charge = 2, TRUE
							   my_peak_mzs,
							   mass_iso_window,
							   2));
	assert(!is_c13_isotopologue(160.50167744, // Outside window charge = 2, FALSE
							   my_peak_mzs,
							   mass_iso_window,
							   2));
// Charge of 2 should not work for 1/3 the window
	assert(!is_c13_isotopologue(160.33445161, // Outside window charge = 2, FALSE
							   my_peak_mzs,
							   mass_iso_window,
							   2));
// Charge of 2 should work when max_charge is 4
	assert(is_c13_isotopologue(160.50167742, // Inside window charge = 2, TRUE
							   my_peak_mzs,
							   mass_iso_window,
							   4));
// 2xC13 should be filtered when max is 2
	assert(is_c13_isotopologue(162.00670968,
							   my_peak_mzs,
							   mass_iso_window,
							   1,
							   2));
// 2xC13 should not be filtered when max is 1
	assert(!is_c13_isotopologue(162.00670968,
							   my_peak_mzs,
							   mass_iso_window,
							   1,
							   1));
}

real[real] select_precursor_ions_topn(
	MSXScan[] scans,
	int n,
	real dew,
	real minimum_intensity,
	real mass_isolation_window,
	real full_scan_time,
	bool filter_c13_isotopologues,
	int max_c13_in_isotopologues = 4)
{
/* Outputs list of [RT:mass] of precursors chosen for fragmentation
 * Arguments:
 *  scans - a MSXScan[] of only MS1 scans of the same polarity
 *  n - number of fragmentations per full scan
 *  dew - time between fragmentation selection of a single precursor ion
 *	minimum_intensity - minimum intensity to consider for fragmentaiton
 *  mass_isolation_window - m/z +/- to be considered the same peak
 *  full_scan_time - time between full scans when fragmenting
 *	filter_c13_isotopologues - whether to filter C13 isotopologues
 *  max_c13_in_isotopologues - the maximum number of C13's in isotopologue peaks
 * Returns:
 *  selected - [RT:mass] of precursor ions chosen
 * 
 */
	real[real] selected; // RT:M/Z  (FINAL list)
	real[real] peaks; // M/Z:intensity (all from current scan)
	real[real] top_peaks; // intensity:M/Z (current scan to add to final)
	real rt_offset = full_scan_time / (n + 1); // for fragmnetation time
	real previous_rt = 0;
	real mass_diff_c13 = 1.00335484;
	for(int scan_number = 0; scan_number<scans.length; ++scan_number) // Iterate over all scans
	{
		real rt = scans[scan_number].retention_time;
		if(scan_number > 0 && rt < (previous_rt + full_scan_time))
		{
			continue;
		}
		peaks = scans[scan_number].peaks;
		foreach(mz; peaks.keys.sort) // Iterate over all peaks in current scan
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
				int[] outside_selected_window;
				for(int i=0; i<selected_rts.length; ++i)
				{
					if(
					   // Check 1: DEW
					   (selected_rts[i] < (rt - dew) || // dew expired - keep
					   selected[selected_rts[i]] > (mz + mass_isolation_window) || // mz > selected mz - keep
					   selected[selected_rts[i]] < (mz - mass_isolation_window)) // mz < selected mz - keep
					   && 
					   // Check 2: Isotopologues
					   (!filter_c13_isotopologues || // If not selected, skip to Check 2
					   !is_c13_isotopologue(mz, //  Not likely to remove 2xC13 isotopologues
						   	 				selected[selected_rts[i]], 
											mass_isolation_window,
						  				    4,
											max_c13_in_isotopologues))
					   )
					{
						outside_selected_window ~= i; //keep
					}
				}
				if(outside_selected_window.length < selected_rts.length) // If not all pass, skip
				{
					outside_selected_window = [];
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
		previous_rt = rt;
	}
	return selected;
}
unittest
{
	import std.stdio;
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
	fourth.retention_time = 9;
	real[real] fourth_peaks = [
		50.0: 1000000.0,
		50.99335485: 1000000.0,
		51.00335484: 1000000.0,
		51.013354829: 1000000.0,
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
		110.0: 1000000.0,
	];
	fourth.peaks = fourth_peaks;
	MSXScan fifth = new MSXScan;
	fifth.level=1;
	fifth.retention_time = 10;
	real[real] fifth_peaks = [
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
	fifth.peaks = fifth_peaks;
	MSXScan[] all_scans = [first, second, third, fourth, fifth];
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
								  full_scan_time,
								  false);
	real[] precursor_rts = my_precursors.keys;
	int num_per_scan = 0;
	foreach(rt; precursor_rts)
	{
		if(rt > 4 && rt < 5)
		{
			++num_per_scan;
		}
	}
	assert(num_per_scan == n); // tests n
	bool same_within_dew = false;
	bool rt_9_mz_110 = false;
	stderr.writeln(my_precursors);	
	foreach(rt, mz; my_precursors)
	{
		if(mz == 110.0 &&
		   rt > 2 &&
		   rt < 5)
		{
			same_within_dew = true;
		}
		if(mz == 110.0 &&
		   rt >= 9)
		{
			stderr.writeln("setting true");
			rt_9_mz_110 = true;
		}
	}
	assert(!same_within_dew); // tests dew
	assert(rt_9_mz_110); // tests exiting dew
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
	real[real] my_other_precursors = select_precursor_ions_topn(all_scans,
										n,
										dew,
										min_intensity,
										mass_isolation_window,
										3,
										false);
	bool scan_before_instrument_ready = false;
	assert(my_other_precursors.length == 9); //tests full scan time
	MSXScan sixth = new MSXScan;
	sixth.level=1;
	sixth.retention_time = 20;
	real[real] sixth_peaks = [
		50.0: 1000000.0,
		50.99335485: 1000000.0,
		51.00335484: 1000000.0,
		51.013354829: 1000000.0,
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
	sixth.peaks = sixth_peaks;
	all_scans = [first, second, third, fourth, fifth, sixth];
	my_precursors = select_precursor_ions_topn(all_scans,
								  n,
								  dew,
								  min_intensity,
								  mass_isolation_window, // 0.01
								  full_scan_time,
								  true);

	foreach(rt, mz; my_precursors)
	{
		assert(mz != 51.00335484); // Check filter
		assert(mz != 50.99335485); // Check filter mass_isolation_window
		assert(mz != 51.01335483); // Check filter mass_isolation_window
	}
}
