import sys
sys.path.append('/switchfinder/')


import SwitchFinder.dotbracket_comparisons as dotbracket_comparisons
import SwitchFinder.glob_vars as glob_vars
import SwitchFinder.utils as utils
import SwitchFinder.check_base_pairing as check_base_pairing
import textwrap
import numpy as np
import copy


class fragment_collection:
    def __init__(self):
        self.body_dict = {}

    def identify_stretches(self,
                           MAX_SPACING = 1,
                           MIN_STEM_LENGHT = 3,
                           MAXIMAL_OVERHANG = 2
                           ):
        for i in self.body_dict:
            self.body_dict[i].identify_stretches(MAX_SPACING, MIN_STEM_LENGHT, MAXIMAL_OVERHANG)
            self.body_dict[i].convert_to_numpy()


    def remove_outer_stretches(self):
        for i in self.body_dict:
            self.body_dict[i].remove_outer_stretches()


    def print(self, num_to_print = 2, width = 91,
              show_changing_structures = True,
              show_unchanged_structures=True,
              show_inner_stems = False,
              show_outer_stems = False,
              do_return = True
              ):
        all_fragments_names = sorted(list(self.body_dict.keys()))
        strings_to_print_list = []
        strings_to_print_list.append("Printing %d first fragments as specified\n" % num_to_print)
        if not do_return:
            print("\n".join(strings_to_print_list))
            strings_to_print_list = []
        for i in range(num_to_print):
            fr_string = self.body_dict[all_fragments_names[i]].print(
                                                         show_changing_structures,
                                                         show_unchanged_structures,
                                                         show_inner_stems,
                                                         show_outer_stems,
                                                         width=width,
                                                         do_return = do_return
                                                         )
            strings_to_print_list.append(fr_string)
        if do_return:
            return "\n".join(strings_to_print_list)


class fragment:
    def __init__(self):
        self.name = ''
        self.major_conf = conformation()
        self.second_conf = conformation()
        self.common_conf = conformation()
        self.sequence = ''
        self.changing_loops_mask = None
        self.A_C_mask = None

    def identify_stretches(self, MAX_SPACING = 1, MIN_STEM_LENGHT = 3, MAXIMAL_OVERHANG = 2):
        self.major_conf.get_all_basepair_stretches()
        self.second_conf.get_all_basepair_stretches()

        self.major_conf.combine_parts_of_the_same_stretch_together(MAX_SPACING)
        self.second_conf.combine_parts_of_the_same_stretch_together(MAX_SPACING)

        unchanged_stretches_major, \
        changed_stretches_major, \
        unchanged_stretches_second, \
        changed_stretches_second = dotbracket_comparisons.remove_changing_parts_both_conformations(
            self.major_conf.stretches_with_spaces, self.second_conf.stretches_with_spaces,
            self.major_conf.string, self.second_conf.string,
            MAXIMAL_OVERHANG = MAXIMAL_OVERHANG
        )
        self.major_conf.changed_stretches = changed_stretches_major
        self.major_conf.unchanged_stretches = unchanged_stretches_major
        self.second_conf.changed_stretches = changed_stretches_second
        self.second_conf.unchanged_stretches = unchanged_stretches_second

        self.major_conf.changed_stretches_no_short = dotbracket_comparisons.remove_short_changing_stems(
            self.major_conf.changed_stretches, MIN_STEM_LENGHT = MIN_STEM_LENGHT
        )
        self.second_conf.changed_stretches_no_short = dotbracket_comparisons.remove_short_changing_stems(
            self.second_conf.changed_stretches, MIN_STEM_LENGHT = MIN_STEM_LENGHT
        )

    def remove_outer_stretches(self):
        no_outer_stretches = dotbracket_comparisons.remove_all_outer_stems(self.major_conf.unchanged_stretches,
                                                    self.major_conf.changed_stretches,
                                                    self.second_conf.changed_stretches)
        self.major_conf.no_outer_stretches = no_outer_stretches
        self.second_conf.no_outer_stretches = no_outer_stretches

        outer_stretches_only = dotbracket_comparisons.two_stretches_difference(self.major_conf.unchanged_stretches,
                                                                             self.major_conf.no_outer_stretches)

        self.major_conf.outer_stretches = outer_stretches_only
        self.second_conf.outer_stretches = outer_stretches_only

        inner_stretches_unchanged = dotbracket_comparisons.two_stretches_difference(self.major_conf.unchanged_stretches,
                                                                                    self.major_conf.outer_stretches)

        self.major_conf.unchanged_inner_stretches = inner_stretches_unchanged
        self.second_conf.unchanged_inner_stretches = inner_stretches_unchanged

    def initialize_pairing_states(self):
        self.major_conf.get_pairing_states()
        self.second_conf.get_pairing_states()


    def get_differential_pairing_states(self):
        if self.major_conf.pairing_states is None or \
                self.second_conf.pairing_states is None:
            self.initialize_pairing_states()

        diff_unpaired_1, diff_unpaired_2 = dotbracket_comparisons.identify_differentially_paired_positions(
            self.major_conf.pairing_states, self.second_conf.pairing_states)
        self.major_conf.differentially_unpaired_states = diff_unpaired_1
        self.second_conf.differentially_unpaired_states = diff_unpaired_2


    def get_constant_pairing_states(self):
        constantly_paired_positions, constantly_unpaired_positions = dotbracket_comparisons.identify_constantly_paired_positions(
            self.major_conf.pairing_states, self.second_conf.pairing_states)
        self.major_conf.constantly_paired_positions = constantly_paired_positions
        self.second_conf.constantly_paired_positions = constantly_paired_positions
        self.major_conf.constantly_unpaired_positions = constantly_unpaired_positions
        self.second_conf.constantly_unpaired_positions = constantly_unpaired_positions


    def get_differential_pairing_states_within_loops(self):
        self.get_differential_pairing_states()
        self.make_mask_of_differential_intervals()
        self.major_conf.differentially_unpaired_states_within_loops = np.logical_and(
                        self.major_conf.differentially_unpaired_states,
                        self.changing_loops_mask)
        self.second_conf.differentially_unpaired_states_within_loops = np.logical_and(
                        self.second_conf.differentially_unpaired_states,
                        self.changing_loops_mask)

    def get_differential_ACs(self):
        self.get_differential_pairing_states_within_loops()
        self.get_AC_mask()
        self.get_constant_pairing_states()
        self.major_conf.differentially_unpaired_ACs = np.logical_and(
                        self.major_conf.differentially_unpaired_states_within_loops,
                        self.A_C_mask)
        self.second_conf.differentially_unpaired_ACs = np.logical_and(
                        self.second_conf.differentially_unpaired_states_within_loops,
                        self.A_C_mask)
        self.major_conf.constant_paired_ACs = np.logical_and(self.major_conf.constantly_paired_positions,
                        self.A_C_mask)
        self.major_conf.constant_unpaired_ACs = np.logical_and(self.major_conf.constantly_unpaired_positions,
                        self.A_C_mask)
        self.second_conf.constant_paired_ACs = self.major_conf.constant_paired_ACs
        self.second_conf.constant_unpaired_ACs = self.major_conf.constant_unpaired_ACs


    def get_differential_intervals(self):
        all_changing_stretches = np.concatenate((
            self.major_conf.changed_stretches_no_short,
            self.second_conf.changed_stretches_no_short),
            axis=0)
        left_halves = all_changing_stretches[:, 0:2]
        right_halves = all_changing_stretches[:, 2:4]
        all_changing_intervals = np.concatenate((left_halves, right_halves), axis=0)
        all_changing_intervals = np.sort(all_changing_intervals, axis=0)
        all_changing_intervals[:, 1] += 1 # the changing stretches coordinates are exclusive, unlike python coordinates
        merged_changing_intervals = utils.merge_intervals(all_changing_intervals)
        return merged_changing_intervals


    def make_mask_of_differential_intervals(self):
        merged_changing_intervals = self.get_differential_intervals()
        changing_state = np.zeros(len(self.sequence), dtype=bool)
        for i in range(merged_changing_intervals.shape[0]):
            curr_interval = merged_changing_intervals[i]
            changing_state[curr_interval[0] : curr_interval[1]] = 1
        self.changing_loops_mask = changing_state


    def get_AC_mask(self):
        A_C_mask = np.zeros(len(self.sequence), dtype=bool)

        for i in range(len(self.sequence)):
            curr_nt = self.sequence[i]
            if curr_nt == 'A' or curr_nt =='C':
                A_C_mask[i] = 1
        self.A_C_mask = A_C_mask


    def convert_to_numpy(self):
        self.major_conf.convert_to_numpy()
        self.second_conf.convert_to_numpy()


    def print(self,
              show_changing_structures = False,
              show_unchanged_structures = False,
              show_inner_stems = True,
              show_outer_stems = True,
              width = 113,
              do_return = False
              ):
        strings_to_print_list = []
        strings_to_print_list.append(self.name)
        strings_to_print_list.append("")
        strings_to_print_list.append("Major loop folding:")
        if not do_return:
            print("\n".join(strings_to_print_list))
            strings_to_print_list = []
        major_str = self.major_conf.print(
                              show_changing_structures,
                              show_unchanged_structures,
                              show_inner_stems,
                              show_outer_stems,
                              width,
                              do_return = do_return
                              )
        strings_to_print_list.append(major_str)
        strings_to_print_list.append("")
        strings_to_print_list.append("Second loop folding:")
        if not do_return:
            print("\n".join(strings_to_print_list))
            strings_to_print_list = []
        second_str = self.second_conf.print(
                               show_changing_structures,
                               show_unchanged_structures,
                               show_inner_stems,
                               show_outer_stems,
                               width,
                               do_return=do_return
                               )
        strings_to_print_list.append(second_str)
        strings_to_print_list.append('\n\n')
        if not do_return:
            print("\n".join(strings_to_print_list))
        else:
            return "\n".join(strings_to_print_list)

    def print_differential_ACs(self):
        major_string = "".join(["|" if x else "." for x in self.major_conf.differentially_unpaired_ACs])
        second_string = "".join(["|" if x else "." for x in self.second_conf.differentially_unpaired_ACs])
        return (major_string, second_string)

    def print_constant_ACs(self):
        paired_string = "".join(["|" if x else "." for x in self.major_conf.constant_paired_ACs])
        unpaired_string = "".join(["|" if x else "." for x in self.second_conf.constant_unpaired_ACs])
        return (paired_string, unpaired_string)


class conformation:
    def __init__(self):
        self.string = ''
        self.pairing_states = None
        self.pairing_states_with_stem_ends = None
        self.numpy = None
        self.constantly_paired_positions = None
        self.constantly_unpaired_positions = None
        self.differentially_unpaired_states = None
        self.differentially_unpaired_states_within_loops = None
        self.differentially_unpaired_ACs = None
        self.constant_paired_ACs = None
        self.constant_unpaired_ACs = None
        self.all_stretches = None
        self.stretches_with_spaces = None
        self.changed_stretches = None
        self.changed_stretches_no_short = None
        self.unchanged_stretches = None
        self.unchanged_inner_stretches = None
        self.outer_stretches = None


    def get_pairing_states(self):
        self.pairing_states = dotbracket_comparisons.convert_string_to_pairing_states(self.string)

    # def get_pairing_states_accounting_for_stem_ends(self):
    #     self.get_pairing_states()
    #     self.pairing_states_with_stem_ends = self.pairing_states.copy()
    #     self.pairing_states_with_stem_ends = dotbracket_comparisons.convert_string_to_pairing_states_with_stem_ends(
    #                     self.pairing_states_with_stem_ends,
    #                     self.all_stretches)


    def convert_to_numpy(self):
        self.numpy = dotbracket_comparisons.convert_string_to_numpy(self.string)


    def get_all_basepair_stretches(self):
        self.all_stretches = dotbracket_comparisons.list_all_stretches_numpy(self.string)


    def combine_parts_of_the_same_stretch_together(self, MAX_SPACING):
        if self.all_stretches is not None:
            self.stretches_with_spaces = dotbracket_comparisons.find_parts_of_the_same_stretch(self.string,
                                                                                               self.all_stretches,
                                                                                               MAX_SPACING = MAX_SPACING)
        else:
            raise Exception("You should find all stretches before finding the stretches with spaces!")


    def print(self,
              show_changing_structures,
              show_unchanged_structures,
              show_inner_stems,
              show_outer_stems,
              width,
              do_return = False
              ):

        strings_to_print_list = []
        strings_to_print_list.append("All stretches")
        strings_to_print_list.append(textwrap.fill(self.string, width=width))


        if show_unchanged_structures:
            if self.unchanged_stretches is not None:
                string_to_print = self.get_unchanged_string()
                strings_to_print_list.append("The stretches that are the same between two conformations")
                strings_to_print_list.append(textwrap.fill(string_to_print, width=width))


        if show_changing_structures:
            if self.changed_stretches_no_short is not None:
                string_to_print = self.get_changing_string()
                strings_to_print_list.append("The stretches that change between two conformations")
                strings_to_print_list.append(textwrap.fill(string_to_print, width=width))

        if show_inner_stems:
            if self.unchanged_inner_stretches is not None:
                string_to_print = self.get_unchanged_inner_stretches()
                strings_to_print_list.append("The INNER stretches that are the same between two conformations")
                strings_to_print_list.append(textwrap.fill(string_to_print, width=width))

        if show_outer_stems:
            if self.outer_stretches is not None:
                string_to_print = self.get_outer_stretches()
                strings_to_print_list.append("The OUTER stretches that are the same between two conformations")
                strings_to_print_list.append(textwrap.fill(string_to_print, width=width))

        final_string_to_print = "\n".join(strings_to_print_list)

        if do_return:
            return final_string_to_print
        else:
            print(final_string_to_print)
            return ""


    def get_changing_string(self):
        if self.changed_stretches_no_short is not None:
            string_to_print = dotbracket_comparisons.string_and_stretches_to_print(
                self.string, self.changed_stretches_no_short
            )
            return string_to_print


    def get_unchanged_string(self):
        if self.changed_stretches_no_short is not None:
            string_to_print = dotbracket_comparisons.string_and_stretches_to_print(
                self.string, self.unchanged_stretches
            )
            return string_to_print


    def get_unchanged_inner_stretches(self):
        if self.unchanged_inner_stretches is not None:
            string_to_print = dotbracket_comparisons.string_and_stretches_to_print(
                self.string, self.unchanged_inner_stretches
            )
            return string_to_print


    def get_outer_stretches(self):
        if self.outer_stretches is not None:
            string_to_print = dotbracket_comparisons.string_and_stretches_to_print(
                self.string, self.outer_stretches
            )
            return string_to_print


class loop:
    def __init__(self):
        self.left_start = None
        self.left_end = None
        self.right_start = None
        self.right_end = None
        self.length = None
        self.left_seq = None
        self.right_seq = None


    def initialize_from_dict(self, inp_dict):
        self.left_start = int(inp_dict['forward'][0])
        self.left_end = int(inp_dict['forward'][1]) + 1
        self.right_start = int(inp_dict['backward'][0])
        self.right_end = int(inp_dict['backward'][1]) + 1
        self.length = self.left_end - self.left_start


    def fill_sequence(self, sequence_np):
        self.left_seq = sequence_np[self.left_start : self.left_end]
        self.right_seq = sequence_np[self.right_start: self.right_end]


    def get_fraction_unpaired(self, allow_canonical = True):
        n_mismatches = check_base_pairing.n_not_paired_bases(
                            self.left_seq, self.right_seq,
                           orientation="reversed",
                           allow_canonical = allow_canonical)
        fraction_unpaired = n_mismatches / self.length
        return fraction_unpaired


    def check_if_complementary(self, allow_canonical = True):
        fraction_unpaired = self.get_fraction_unpaired(allow_canonical = allow_canonical)
        is_complementary = fraction_unpaired == 0
        return is_complementary


    def introduce_random_sequence(self, seed = None,
                                  GC_only = False):
        if not seed is None:
            np.random.seed(seed)
        if not GC_only:
            low = 1
            high = 5
        else:
            low = 2
            high = 4
        random_sequence_array = np.random.randint(low=low, high=high, size=self.length)
        self.left_seq = random_sequence_array
        self.right_seq = check_base_pairing.make_reverse_complement(self.left_seq)


    def change_one_side(self, change_left_side = True, seed = None):
        if not seed is None:
            np.random.seed(seed)
        random_sequence_array = np.random.randint(low=1, high=5, size=self.length)
        if change_left_side:
            self.left_seq = random_sequence_array
        else:
            self.right_seq = random_sequence_array


    def change_one_side_given_seq(self,
                                  new_seq_array,
                                  left_subset_border,
                                  right_subset_border,
                                  change_left_side = True,
                                  make_complementary = False):
        subset_length = right_subset_border - left_subset_border
        compl_part_beginning = self.length - right_subset_border
        compl_part_end = compl_part_beginning + subset_length

        if change_left_side:
            self.left_seq = new_seq_array
            if make_complementary:
                subset_to_make_complementary = self.left_seq[left_subset_border : right_subset_border]
                complementary_part = check_base_pairing.make_reverse_complement(subset_to_make_complementary)
                self.right_seq[compl_part_beginning : compl_part_end] = complementary_part
        else:
            self.right_seq = new_seq_array
            if make_complementary:
                subset_to_make_complementary = self.right_seq[left_subset_border : right_subset_border]
                complementary_part = check_base_pairing.make_reverse_complement(subset_to_make_complementary)
                self.left_seq[compl_part_beginning : compl_part_end] = complementary_part

    def mask_out(self):
        mask_none_array = np.full(self.length, glob_vars._none)
        self.left_seq = mask_none_array
        self.right_seq = mask_none_array


    def longest_paired_stretch(self,
                               allow_canonical = True):
        max_bp_stretch = check_base_pairing.max_base_pairs_in_a_row(
                            self.left_seq, self.right_seq,
                            allow_canonical=allow_canonical)
        return max_bp_stretch


    def change_original_sequence_accordingly(self, sequence_np):
        out_sequence = sequence_np.copy()
        out_sequence[self.left_start : self.left_end] = self.left_seq
        out_sequence[self.right_start: self.right_end] = self.right_seq
        return out_sequence


    def scan_original_sequence_for_bp_stretches(self, sequence_np, allow_canonical = False):
        # scans the provided sequence with both sides of the stem and checks for the longest possible base pairing
        # masks out the stem itself
        masked_sequence_np = copy.deepcopy(sequence_np)
        masked_sequence_np[self.left_start : self.left_end] = -1
        masked_sequence_np[self.right_start: self.right_end] = -1
        max_stretch = self.scan_other_sequence_for_bp_stretches(masked_sequence_np,
                                                  allow_canonical = allow_canonical)
        return max_stretch


    def scan_original_sequence_for_bp_stretches_one_loop(self, sequence_np, allow_canonical = False,
                                                      check_left_side = True):
        masked_sequence_np = copy.deepcopy(sequence_np)
        masked_sequence_np[self.left_start : self.left_end] = -1
        masked_sequence_np[self.right_start: self.right_end] = -1
        max_stretch = self.scan_other_sequence_for_bp_stretches_one_loop(
                                                        sequence_np,
                                                        allow_canonical = allow_canonical,
                                                        check_left_side = check_left_side)
        return max_stretch


    def scan_other_sequence_for_bp_stretches(self, sequence_np, allow_canonical = False):
        # scans the provided sequence with both sides of the stem and checks for the longest possible base pairing
        # array of 2 x N holding lengths of the maximal stretches in each possible pairing
        # the first dimension is of size 2 since I am scanning both forward and reverse sequence

        max_stretch_left = self.scan_other_sequence_for_bp_stretches_one_loop(
                                                        sequence_np,
                                                        allow_canonical = allow_canonical,
                                                        check_left_side = True)
        max_stretch_right = self.scan_other_sequence_for_bp_stretches_one_loop(
                                                        sequence_np,
                                                        allow_canonical = allow_canonical,
                                                        check_left_side = False)
        max_stretch = max(max_stretch_left, max_stretch_right)
        return max_stretch


    def scan_other_sequence_for_bp_stretches_one_loop(self, sequence_np, allow_canonical = False,
                                                      check_left_side = True):
        # scans the provided sequence with one side of the stem and checks for the longest possible base pairing
        # array of 1 x N holding lengths of the maximal stretches in each possible pairing
        full_seq_max_bp_stretches = np.zeros(sequence_np.shape[0] - self.length + 1)

        if check_left_side:
            query_sequence = self.left_seq
        else:
            query_sequence = self.right_seq

        for i in range(full_seq_max_bp_stretches.shape[0]):
            full_seq_max_bp_stretches[i] = check_base_pairing.max_base_pairs_in_a_row(
                                query_sequence,
                                sequence_np[i: i + self.length],
                                allow_canonical = allow_canonical)
        max_stretch = np.max(full_seq_max_bp_stretches)
        return max_stretch


    def scan_original_sequence_for_number_of_bps(self, sequence_np, allow_canonical = False):
        # scans the provided sequence with both sides of the stem and checks for the maximal possible number of base pairs
        # masks out the stem itself
        masked_sequence_np = copy.deepcopy(sequence_np)
        masked_sequence_np[self.left_start : self.left_end] = -1
        masked_sequence_np[self.right_start: self.right_end] = -1
        max_number_of_base_pairs = self.scan_other_sequence_for_number_of_bps(masked_sequence_np,
                                                  allow_canonical = allow_canonical)
        return max_number_of_base_pairs


    def scan_original_sequence_for_number_of_bps_one_loop(self, sequence_np, allow_canonical = False,
                                                      check_left_side = True):
        masked_sequence_np = copy.deepcopy(sequence_np)
        masked_sequence_np[self.left_start : self.left_end] = -1
        masked_sequence_np[self.right_start: self.right_end] = -1
        max_number_of_base_pairs = self.scan_other_sequence_for_number_of_bp_one_loop(
                                                        sequence_np,
                                                        allow_canonical = allow_canonical,
                                                        check_left_side = check_left_side)
        return max_number_of_base_pairs


    def scan_other_sequence_for_number_of_bps(self, sequence_np, allow_canonical = False):
        # scans the provided sequence with both sides of the stem and checks for the maximal possible number of base pairs
        # array of 2 x N holding lengths of the maximal possible number of base pairs in each position
        # the first dimension is of size 2 since I am scanning both forward and reverse sequence

        max_number_of_base_pairs_left = self.scan_other_sequence_for_number_of_bp_one_loop(
                                                        sequence_np,
                                                        allow_canonical = allow_canonical,
                                                        check_left_side = True)
        max_number_of_base_pairs_right = self.scan_other_sequence_for_number_of_bp_one_loop(
                                                        sequence_np,
                                                        allow_canonical = allow_canonical,
                                                        check_left_side = False)
        max_number_of_base_pairs = max(max_number_of_base_pairs_left, max_number_of_base_pairs_right)
        return max_number_of_base_pairs


    def scan_other_sequence_for_number_of_bp_one_loop(self, sequence_np, allow_canonical = False,
                                                      check_left_side = True):
        # scans the provided sequence with one side of the stem and counts the number of base pairs
        # array of 1 x N holding number of possible base pairs at each positions
        number_of_base_pairs = np.zeros(sequence_np.shape[0] - self.length + 1)
        if check_left_side:
            query_sequence = self.left_seq
        else:
            query_sequence = self.right_seq

        for i in range(number_of_base_pairs.shape[0]):
            number_of_mismatches = check_base_pairing.n_not_paired_bases(
                                                query_sequence,
                                                sequence_np[i: i + self.length],
                                                allow_canonical=allow_canonical)
            number_of_base_pairs[i] = self.length - number_of_mismatches
        max_number_of_base_pairs = np.max(number_of_base_pairs)
        return max_number_of_base_pairs


    def compare_to_base_pair_array(self, base_pair_obj,
                                   do_print = False):
        bps_of_interest = base_pair_obj.array[self.left_start : self.left_end]
        expected_bps = np.arange(self.right_end - 1, self.right_start - 1, -1)
        same_bps_number = (expected_bps == bps_of_interest).sum()
        if do_print:
            print("bps of interest")
            print(bps_of_interest)
            print("expected bps")
            print(expected_bps)

        return same_bps_number


    def copy(self):
        new_loop = copy.deepcopy(self)
        return new_loop


    def print(self, do_return = False):
        string_to_print = "Left stretch: %d - %d ; right stretch: %d - %d\n" % \
                (self.left_start, self.left_end, self.right_start, self.right_end)
        if not self.left_seq is None and not self.right_seq is None:
            string_to_print += "Left sequence: %s ; right sequence %s\n" % \
               (utils.array_to_string(self.left_seq), utils.array_to_string(self.right_seq))
        if do_return:
            return string_to_print
        print(string_to_print)


    def get_base_pairing_array(self, sequence_np):
        base_pairing_state_array = np.full(sequence_np.shape[0], glob_vars._loop, dtype = np.uint8)
        base_pairing_state_array[self.left_start : self.left_end] = glob_vars._left_stem
        base_pairing_state_array[self.right_start: self.right_end] = glob_vars._right_stem
        return base_pairing_state_array


    def get_dot_bracket_string(self, sequence_np):
        base_pairing_state_array = self.get_base_pairing_array(sequence_np)
        base_pairing_state_strings_array = [glob_vars._extended_structure_to_char[x] for x in base_pairing_state_array]
        base_pairing_state_string = "".join(base_pairing_state_strings_array)
        return base_pairing_state_string



class base_pair_array:
    def __init__(self, inp_string):
        self.string = inp_string
        self.fill_array_from_string()

    def fill_array_from_string(self, unpaired_value = -1):
        bp_array = np.full(len(self.string), unpaired_value, dtype = np.int)
        stretches_np = dotbracket_comparisons.list_all_stretches_numpy(self.string)
        for i in range(stretches_np.shape[0]):
            current_stretch = stretches_np[i, :]
            stretch_length = current_stretch[1] - current_stretch[0] + 1
            for k in range(stretch_length):
                left_position = current_stretch[0] + k
                right_position = current_stretch[3] - k
                bp_array[left_position] = right_position
                bp_array[right_position] = left_position
        self.array = bp_array

    def array_to_string(self):
        pairing_mode_array = np.zeros(self.array.shape[0], dtype = np.uint8)
        indices = np.arange(self.array.shape[0])
        pairing_mode_array[indices < self.array] = glob_vars._left_stem
        pairing_mode_array[indices > self.array] = glob_vars._right_stem
        pairing_mode_array[self.array < 0] = glob_vars._loop
        string_to_write_list = [glob_vars._extended_structure_to_char[x] for x in pairing_mode_array]
        string_to_write = "".join(string_to_write_list)
        return string_to_write





