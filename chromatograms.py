# classes for peak chromatograms
import bisect
import numpy as np
import pylab as plt
from scipy.stats import pearsonr

class PeakChromatogram(object):
    def __init__(self):
        self.mz_list = []
        self.rt_list = []
        self.intensity_list = []
        self.scan_n_list = []
        self.n_points = len(self.mz_list)

    def mean_mz(self):
        return sum(self.mz_list)/self.n_points

    def add_point(self,mz,rt,intensity,scan_number):
        self.mz_list.append(mz)
        self.rt_list.append(rt)
        self.intensity_list.append(intensity)
        self.scan_n_list.append(scan_number)
        self.n_points += 1

    def plot(self,**kwargs):
        plt.figure(**kwargs)
        plt.plot(self.rt_list,self.intensity_list)

    def correlation(self,other):
        if self.n_points < 2 or other.n_points < 2:
            return 0.0 
        if other.scan_n_list[0] > self.scan_n_list[-1] or other.scan_n_list[-1] < self.scan_n_list[0]:
            return 0.0 # no overlap
        min_scan = min(self.scan_n_list[0],other.scan_n_list[0])
        max_scan = max(self.scan_n_list[-1],other.scan_n_list[-1])
        scans = range(min_scan,max_scan+1)
        self_intensity = []
        other_intensity = []
        for s in scans:
            try:
                self_pos = self.scan_n_list.index(s)
                self_intensity.append(self.intensity_list[self_pos])
            except:
                self_intensity.append(0.0)
            try:
                other_pos = other.scan_n_list.index(s)
                other_intensity.append(other.intensity_list[other_pos])
            except:
                other_intensity.append(0.0)

        c,p = pearsonr(self_intensity,other_intensity)
        return c            
        

def add_chromatograms_to_boxes(boxes,mz_file_obj):
    for box in boxes:
        box.peak_chromatogram = PeakChromatogram()
    
    boxes.sort(key = lambda x: x.rt_range[0])
    remaining_boxes = [b for b in boxes]
    completed_boxes = []
    current_boxes = []
    for scan in mz_file_obj.scans:
        if scan.ms_level > 1:
            continue
        scan_rt_seconds = scan.rt_in_seconds
        # move any remaining boxes that are now live into current
        new_boxes = list(filter(lambda x: scan_rt_seconds >= x.rt_range_in_seconds[0],remaining_boxes))
        current_boxes += new_boxes
        for box in new_boxes:
            pos = remaining_boxes.index(box)
            del remaining_boxes[pos]
        # remove any current that are no longer
        done_boxes = list(filter(lambda x: scan_rt_seconds > x.rt_range_in_seconds[1],current_boxes))
        completed_boxes += done_boxes
        for box in done_boxes:
            pos = current_boxes.index(box)
            del current_boxes[pos]
        
        if len(current_boxes) > 0:
            mz,intensity = zip(*scan.peaks)
            # find peaks for each box in this scan (use the most intense if > 1)
            for box in current_boxes:
                assert box.rt_range_in_seconds[0] <= scan_rt_seconds
                assert box.rt_range_in_seconds[1] >= scan_rt_seconds
                min_mz,max_mz = box.mz_range
                lpos = bisect.bisect_left(mz,min_mz)
                rpos = bisect.bisect_right(mz,max_mz)
                if lpos == rpos:
                    continue # no peak
                mpos = np.argmax(intensity[lpos:rpos])
                mpos += lpos
                box.peak_chromatogram.add_point(mz[mpos],scan_rt_seconds,intensity[mpos],scan.scan_no)
                # print(len(box.peak_chromatogram.mz_list))

        # print(len(completed_boxes),len(current_boxes),len(remaining_boxes))

def cluster_box_chromatograms(boxes,threshold = 0.7,max_mz_diff = 100):
    correlations = {box:{} for box in boxes}
    for i,box in enumerate(boxes[:-1]):
        for box2 in boxes[i+1:]:
            c = box.peak_chromatogram.correlation(box2.peak_chromatogram)
            if c > 0:
                correlations[box][box2] = c
                correlations[box2][box] = c
    
    groups = []
    remaining_boxes = [b for b in boxes]

    while len(remaining_boxes) > 0:
        remaining_boxes = list(remaining_boxes)
        remaining_boxes.sort(key = lambda x: x.height,reverse = True)
        new_group = set()
        seed_box = remaining_boxes[0]
        new_group.add(seed_box)
        remaining_boxes = set(remaining_boxes)
        remaining_boxes.remove(seed_box)
        for box,c in correlations[seed_box].items():
            if box in remaining_boxes and c >= threshold and abs(seed_box.mz - box.mz) <= max_mz_diff:
                new_group.add(box)
                remaining_boxes.remove(box)
        groups.append(new_group)


    return groups,correlations

def peak_group(boxes,seed_box,threshold = 0.7,max_mz_diff = 100):
    # find all the boxes that link to a single box, with a higher mz
    # all the ones with a greater mz
    sub_boxes = list(filter(lambda x: x.mz > seed_box.mz and x.mz < seed_box.mz + max_mz_diff,boxes))
    # compute correlations
    correlations = []
    for box in sub_boxes:
        co = seed_box.peak_chromatogram.correlation(box.peak_chromatogram)
        if co > 0:
            correlations.append((box,co))
    
    correlations = set(filter(lambda x: x[1] >= threshold,correlations))
    correlations.add(seed_box)

    return correlations
    