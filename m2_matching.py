# code for ms2 processing
import pymzml
class MZMLScan(object):
    def __init__(self,scan_no,source_file,ms_level,peaks,rt_in_minutes,precursor_mz = None):
        self.scan_no = scan_no
        self.source_file = source_file
        self.ms_level = ms_level
        self.peaks = peaks
        self.rt_in_minutes = rt_in_minutes
        self.precursor_mz = self.precursor_mz

class MZMLFile(object):
    def __init__(self,file_name):
        self.file_name = file_name
        self.scans = []
        self._load_file()
    
    def _load_file(self):
        run = pymzml.reader.Run(self.file_name)
        scan_no = 0
        for scan in run:
            rt = scan.scan_time_in_minutes()
            ms_level = scan.ms_level
            peaks = scan.peaks('centroided')
            if ms_level == 2:
                precursor_mz = scan.selected_precursors[0]['mz']
            else:
                precursor_mz = None
            scan_no += 1
            self.scans.append(MZMLScan(scan_no,self.file_name,
                              ms_level,peaks,rt,precursor_mz))
