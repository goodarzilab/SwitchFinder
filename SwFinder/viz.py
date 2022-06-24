import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
import scipy.stats
from . import glob_vars
from . import perturbations


def check_the_correlations(inp_df, colnames_1, colnames_2, n_cols=2):
    counter = 0
    fig, ax = plt.subplots(nrows=len(colnames_1) // n_cols, ncols=n_cols)
    fig.set_size_inches(8, 16)

    for c1, c2 in zip(colnames_1, colnames_2):
        bin_name = c1.split('_')[1]
        assert (bin_name == c2.split('_')[1])
        counter += 1
        current_correlation, pv = scipy.stats.spearmanr(inp_df[c1], inp_df[c2])
        title_plot = "Bin %s; spearman corr %.2f" % (bin_name, round(current_correlation, 2))
        plt.subplot(len(colnames_1) // n_cols, n_cols, counter)
        plt.scatter(inp_df[c1], inp_df[c2], s=0.1, alpha=0.3)
        plt.title(title_plot)
        plt.xlim(0, np.quantile(inp_df[c1], 0.99))
        plt.ylim(0, np.quantile(inp_df[c2], 0.99))
        plt.tick_params(
            axis='both',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False,
            labeltop=False,
            left=False,  # ticks along the bottom edge are off
            right=False,  # ticks along the top edge are off
            labelleft=False,
            labelright=False
        )
        plt.tick_params(
            axis='y',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False,
            labeltop=False
        )

    plt.show()


def scale_shape_values_to_RGB(shape_array, colormap):
    out_colors = np.zeros((shape_array.shape[0], 3), dtype=np.float)
    vmin = glob_vars._SHAPE_VMIN
    vmax = glob_vars._SHAPE_VMAX
    normalizer = Normalize(vmin, vmax)
    norm_colormap = cm.ScalarMappable(norm=normalizer, cmap=colormap)

    new_array = shape_array.copy()
    sub_array = new_array[shape_array >= 0]
    np.clip(sub_array, vmin, vmax)

    for i in range(shape_array.shape[0]):
        curr_value = shape_array[i]
        if curr_value < vmin:
            continue
        if curr_value > vmax:
            curr_value = vmax
        # from here https://stackoverflow.com/questions/15140072/how-to-map-number-to-color-using-matplotlibs-colormap
        r, g, b, alpha = norm_colormap.to_rgba(curr_value)
        out_colors[i, 0] = r
        out_colors[i, 1] = g
        out_colors[i, 2] = b
    return out_colors


# https://stackoverflow.com/questions/29779079/adding-a-scatter-of-points-to-a-boxplot-using-matplotlib
def dms_values_boxplot(array_1, array_2, ax = None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    sns.violinplot(data = [array_1, array_2], ax = ax, color=".8", **kwargs)
    sns.stripplot(data = [array_1, array_2], ax = ax, **kwargs)
    return ax

def two_correlations_scatterplot(ax, x, y,
                                 xlim = (-1, 1),
                                 ylim = (-1, 1)):
    ax.scatter(x, y, s=5)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.grid()
    y_gridlines = ax.yaxis.get_gridlines()
    y_gridlines[4].set_color("k")
    y_gridlines[4].set_linewidth(1)
    x_gridlines = ax.xaxis.get_gridlines()
    x_gridlines[4].set_color("k")
    x_gridlines[4].set_linewidth(1)
    return ax

def plot_post_normalization_histograms(inp_df,
                                       column_names,
                                       figsize = (16,20),
                                       range = (0, 2000),
                                       bins = 15):
    fig, ax = plt.subplots(ncols = 2, nrows = 8, figsize = figsize)
    for i, sample in enumerate(column_names):
        ax[i % 8, i // 8].hist(inp_df[sample],
                            bins = bins,
                            range = range)
    plt.show()


def mut_corr_scatter_plots(mut_df,
                           title = '',
                           figsize = (23,8),
                           xlim = (-1, 1),
                           ylim = (-1, 1)):
    fig, axs = plt.subplots(ncols = 3, figsize = figsize)
    axs[0] = two_correlations_scatterplot(axs[0],
                    mut_df['similar_avg_1'],
                    mut_df['dissimilar_avg_1'],
                    xlim = xlim,
                    ylim = ylim)
    axs[0].set_title("%s: correlations rep 1" % (title))
    axs[0].set_xlabel("similar mutations")
    axs[0].set_ylabel("dissimilar mutations")
    axs[1] = two_correlations_scatterplot(axs[1],
                    mut_df['similar_avg_2'],
                    mut_df['dissimilar_avg_2'],
                    xlim = xlim,
                    ylim = ylim)
    axs[1].set_title("%s: correlations rep 2" % (title))
    axs[1].set_xlabel("similar mutations")
    axs[1].set_ylabel("dissimilar mutations")
    axs[2] = two_correlations_scatterplot(axs[2],
                    mut_df['similar_avg_both'],
                    mut_df['dissimilar_avg_both'],
                    xlim = xlim,
                    ylim = ylim)
    axs[2].set_title("%s: correlations average replicates" % (title))
    axs[2].set_xlabel("similar mutations")
    axs[2].set_ylabel("dissimilar mutations")
    plt.show()


def scatter_avg_differences(avg_diff_series_1,
                            avg_diff_series_2,
                            ax,
                            title = '',
                            label_1 = '',
                            label_2 = '',
                            xlim=(-1, 1),
                            ylim=(-1, 1)
                            ):
    ax.scatter(avg_diff_series_1,
               avg_diff_series_2,
               s=5)
    ax.set_title(title)
    ax.set_xlabel(label_1)
    ax.set_ylabel(label_2)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.grid()

    return ax


def visualize_one_subarray(ax,
                           array,
                           n_elements,
                           yticks_labels,
                           label,
                           title_fontsize):
    ax.imshow(array[0:n_elements, ],
                  cmap='YlOrRd')
    ax.set_yticks(np.arange(n_elements))
    if not yticks_labels is None:
        ax.set_yticklabels(yticks_labels, fontdict=None, minor=False)
    ax.set_title("Replicate 1: \n %s" % label, fontsize = title_fontsize)
    return ax


def visualize_two_mpra_arrays(array_1_raw,
                             array_2_raw,
                             n_elements,
                             x_dimension=10,
                             do_scale=False,
                             do_return=True,
                             label = "",
                             title_fontsize = 10):
    if do_scale:
        array_1 = array_1_raw / array_1_raw.sum(axis=1)[:, np.newaxis]
        array_2 = array_2_raw / array_2_raw.sum(axis=1)[:, np.newaxis]
    else:
        array_1 = array_1_raw.copy()
        array_2 = array_2_raw.copy()

    y_dimension = x_dimension * n_elements / (2 * array_1.shape[1] + 2)
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(x_dimension, y_dimension))

    axs[0] = visualize_one_subarray(axs[0],
                           array_1,
                           n_elements,
                           label,
                           title_fontsize)
    axs[1] = visualize_one_subarray(axs[1],
                           array_2,
                           n_elements,
                           label,
                           title_fontsize)

    if do_return:
        return fig
    else:
        plt.show()


def visualize_four_mpra_arrays(array_1_raw,
                             array_2_raw,
                             array_3_raw,
                             array_4_raw,
                             n_elements,
                             ytick_labels,
                             x_dimension=10,
                             do_scale=False,
                             do_return=True,
                             label = "",
                             title_fontsize = 10):
    if do_scale:
        array_1 = array_1_raw / array_1_raw.sum(axis=1)[:, np.newaxis]
        array_2 = array_2_raw / array_2_raw.sum(axis=1)[:, np.newaxis]
        array_3 = array_3_raw / array_3_raw.sum(axis=1)[:, np.newaxis]
        array_4 = array_4_raw / array_4_raw.sum(axis=1)[:, np.newaxis]
        # array_1 = (array_1_raw - array_1_raw.min()) / (array_1_raw.max() - array_1_raw.min())
        # array_2 = (array_2_raw - array_2_raw.min()) / (array_2_raw.max() - array_2_raw.min())
        # array_3 = (array_3_raw - array_3_raw.min()) / (array_3_raw.max() - array_3_raw.min())
        # array_4 = (array_4_raw - array_4_raw.min()) / (array_4_raw.max() - array_4_raw.min())
    else:
        array_1 = array_1_raw.copy()
        array_2 = array_2_raw.copy()
        array_3 = array_3_raw.copy()
        array_4 = array_4_raw.copy()

    y_dimension = x_dimension * n_elements / (4 * array_1.shape[1] + 4)
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(x_dimension, y_dimension))

    axs[0] = visualize_one_subarray(axs[0],
                           array_1,
                           n_elements,
                           ytick_labels,
                           label = "gDNA: \n %s" % label,
                           title_fontsize = title_fontsize)
    axs[1] = visualize_one_subarray(axs[1],
                           array_2,
                           n_elements,
                           ytick_labels,
                           label = "gDNA: \n %s" % label,
                           title_fontsize = title_fontsize)
    axs[2] = visualize_one_subarray(axs[2],
                           array_3,
                           n_elements,
                           ytick_labels,
                           label = "RNA: \n %s" % label,
                           title_fontsize = title_fontsize)
    axs[3] = visualize_one_subarray(axs[3],
                           array_4,
                           n_elements,
                           ytick_labels,
                           label = "RNA: \n %s" % label,
                           title_fontsize = title_fontsize)

    if do_return:
        return fig
    else:
        plt.show()


def visualize_pert_profiles_both_replicates(index_selection,
                                            corr_df,
                                            counts_df,
                                            perturbations_order,
                                            replicate_1_column_names,
                                            replicate_2_column_names,
                                            title_fontsize,
                                            do_print = False
                                            ):
    for index, row in corr_df.loc[index_selection].iterrows():
        counts_rep_1, counts_rep_2 = perturbations.subset_counts_five_perturbations(
            row,
            counts_df,
            perturbations_order,
            replicate_1_column_names,
            replicate_2_column_names)
        visualize_two_mpra_arrays(counts_rep_1.to_numpy(),
                                 counts_rep_2.to_numpy(),
                                 counts_rep_1.shape[0],
                                 x_dimension=10,
                                 do_scale=False,
                                 do_return=True,
                                 label = "%s\n%s" %
                                         (corr_df.loc[index, 'gene_name'],
                                          corr_df.loc[index, 'fr_name']),
                                 title_fontsize = title_fontsize)
    if do_print:
        sub_corr_df = corr_df.loc[index_selection]
        for index, row in sub_corr_df.iterrows():
            current_string = "%s\t%s" % (row['gene_name'], row['fr_name'])
            print(current_string)


def visualize_pert_profiles_gDNA_and_RNA(
                                    index_selection,
                                    gdna_corr_df,
                                    gdna_counts_df,
                                    rna_corr_df,
                                    rna_counts_df,
                                    perturbations_order,
                                    replicate_1_column_names,
                                    replicate_2_column_names,
                                    title_fontsize,
                                    x_dimension = 10,
                                    do_print = False,
                                    do_scale = False
                                    ):
    for index in index_selection:
        rna_counts_rep_1, rna_counts_rep_2 = perturbations.subset_counts_five_perturbations(
                        rna_corr_df.loc[index],
                        rna_counts_df,
                        perturbations_order,
                        replicate_1_column_names,
                        replicate_2_column_names)
        gdna_counts_rep_1, gdna_counts_rep_2 = perturbations.subset_counts_five_perturbations(
                        gdna_corr_df.loc[index],
                        gdna_counts_df,
                        perturbations_order,
                        replicate_1_column_names,
                        replicate_2_column_names)
        visualize_four_mpra_arrays(
                                 gdna_counts_rep_1.to_numpy(),
                                 gdna_counts_rep_2.to_numpy(),
                                 rna_counts_rep_1.to_numpy(),
                                 rna_counts_rep_2.to_numpy(),
                                 gdna_counts_rep_1.shape[0],
                                 ytick_labels = None,
                                 x_dimension = x_dimension,
                                 do_scale=do_scale,
                                 do_return=True,
                                 label = "%s\n%s" %
                                         (rna_corr_df.loc[index, 'gene_name'],
                                          rna_corr_df.loc[index, 'fr_name']),
                                 title_fontsize = title_fontsize)
    if do_print:
        sub_corr_df = rna_corr_df.loc[index_selection]
        for index, row in sub_corr_df.iterrows():
            current_string = "%s\t%s\t%s" % (row['gene_name'], row['fr_name'], row['source'])
            print(current_string)

def visualize_flow_perturbation_plots(fragments_list,
                                    perturbations_order,
                                    perturbations_colors_dict,
                                    pert_to_well_dict,
                                    values_dict,
                                    bins,
                                    ncols = 3,
                                    nrows = 7,
                                    figsize=(20, 35),
                                    min_cells = 100,
                                    xlim=(0, 6)):
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=figsize)

    for i, fr in enumerate(fragments_list):
        col = i % ncols
        row = i // ncols
        for pert_name in perturbations_order:
            current_sample_name = "%s_%s" % (fr, pert_name)
            current_well = pert_to_well_dict[current_sample_name]
            current_values = values_dict[current_well]

            # if there aren't even 100 cells
            if len(current_values) < min_cells:
                continue

            curr_distplot = sns.distplot(current_values,
                                         bins=bins,
                                         hist=True,
                                         kde=True,
                                         color=perturbations_colors_dict[pert_name],
                                         ax=ax[row, col])
            curr_distplot.set(xlim = xlim)
            curr_distplot.set_title(fr)


def visualize_bp_pair_matrices_difference_single_mutation(x):
    fig, axs = plt.subplots(nrows = 1, ncols = 3, figsize = (20, 10))
    confl_loops_difference = x['matrix_constrained_second'] - x['matrix_constrained_major']
    axs[0].imshow(confl_loops_difference, cmap = 'bwr',
                  vmin = -np.max(confl_loops_difference),
                  vmax = np.max(confl_loops_difference))
    axs[0].set_title("Conflicting loops difference (%s)" % x['gene_name'])
    mutated_difference = x['matrix_original'] - x['matrix_mutated']
    axs[1].imshow(mutated_difference, cmap = 'bwr',
                  vmin = -np.max(mutated_difference),
                  vmax = np.max(mutated_difference))
    axs[1].set_title("Mutated vs original difference (mutation %s)" % (x['mut_name']))
    rescued_difference = x['matrix_original'] - x['matrix_rescued']
    axs[2].imshow(rescued_difference, cmap = 'bwr',
                  vmin = -np.max(rescued_difference),
                  vmax = np.max(rescued_difference))
    axs[2].set_title("Mutated vs original difference")
    #plt.suptitle("%s; the mutation difference is %.5f" % (title_prefix, rescued_difference.sum()))
    plt.show()


def visualize_bp_pair_matrices_single_mutation(x):
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20, 10))
    vmin = 0
    vmax = np.max(x['matrix_original'])

    axs[0].imshow(x['matrix_original'], cmap='Greys', vmin=vmin, vmax=vmax)
    axs[0].set_title("%s %s %s" % (x['gene_name'], x['fr_name'], x['source']))
    axs[1].imshow(x['matrix_mutated'], cmap='Greys', vmin=vmin, vmax=vmax)
    axs[1].set_title("Mutated (mutation %s)" % (x['mut_name']))
    axs[2].imshow(x['matrix_rescued'], cmap='Greys', vmin=vmin, vmax=vmax)
    axs[2].set_title("Rescued")
    plt.show()


def visualize_mutation_prob_ratios(mut_prob_ratios_df,
                                   name,
                                   figsize = (15, 7),
                                   s = 5,
                                   alpha = 0.3):
    fig, axs = plt.subplots(ncols = 2, figsize = figsize)
    axs[0].scatter(mut_prob_ratios_df['mutated_first_loop_ratio'],
                   mut_prob_ratios_df['mutated_second_loop_ratio'],
                   s = s, alpha = alpha)
    axs[0].set_title("%s - mutation" % (name))
    axs[0].set_xlabel('first loop ratio')
    axs[0].set_ylabel('second loop ratio')
    axs[1].scatter(mut_prob_ratios_df['rescued_first_loop_ratio'],
                   mut_prob_ratios_df['rescued_second_loop_ratio'],
                   s = s, alpha = alpha)
    axs[1].set_title("%s - rescued" % (name))
    axs[1].set_xlabel('first loop ratio')
    axs[1].set_ylabel('second loop ratio')
    plt.show()


def visualize_mutation_prob(mut_prob_ratios_df,
                                   name,
                                   figsize = (15, 7),
                                   s = 5,
                                   alpha = 0.3):
    fig, axs = plt.subplots(ncols = 2, figsize = figsize)
    axs[0].scatter(mut_prob_ratios_df['mutated_first_loop'],
                   mut_prob_ratios_df['mutated_second_loop'],
                   s = s, alpha = alpha)
    axs[0].scatter(mut_prob_ratios_df['original_first_loop'],
                   mut_prob_ratios_df['original_second_loop'],
                   s = s * 1.1, c = "red")
    axs[0].set_title("%s - mutation" % (name))
    axs[0].set_xlabel('first loop probability')
    axs[0].set_ylabel('second loop probability')
    axs[0].set_xlim((-0.03, 1.03))
    axs[0].set_ylim((-0.03, 1.03))
    axs[1].scatter(mut_prob_ratios_df['rescued_first_loop'],
                   mut_prob_ratios_df['rescued_second_loop'],
                   s = s, alpha = alpha)
    axs[1].scatter(mut_prob_ratios_df['original_first_loop'],
                   mut_prob_ratios_df['original_second_loop'],
                   s = s * 1.1, c = "red")
    axs[1].set_title("%s - rescued" % (name))
    axs[1].set_xlabel('first loop probability')
    axs[1].set_ylabel('second loop probability')
    axs[1].set_xlim((-0.03, 1.03))
    axs[1].set_ylim((-0.03, 1.03))
    plt.show()


def visualize_single_mutation_rescue(samples_organization_dict,
                                    title,
                                    conf_order,
                                    conf_colors_dict,
                                    values_dict,
                                    bins,
                                    ax,
                                    xlim=(0, 6),
                                    ):
    for conf_name in conf_order:
        if conf_name not in samples_organization_dict:
            continue
        current_sample_name = samples_organization_dict[conf_name]
        current_values = values_dict[current_sample_name]
        median_value = np.median(current_values)

        curr_distplot = sns.distplot(current_values,
                                     label = current_sample_name,
                                     bins=bins,
                                     hist=True,
                                     kde=True,
                                     color = conf_colors_dict[conf_name],
                                     hist_kws={"range": [xlim[0], xlim[1]]},
                                     ax = ax)
        # add median line
        ax.axvline(median_value, color = conf_colors_dict[conf_name], linestyle='dashed', linewidth=1)
    curr_distplot.set(xlim = xlim)
    curr_distplot.set_title(title, fontsize='large')
    curr_distplot.legend(loc='upper right')



def visualize_pert_profiles_gDNA_and_RNA_v2(
                                    index_selection,
                                    mutation_summary_df,
                                    gdna_counts_df,
                                    rna_counts_df,
                                    perturbations_order,
                                    perturbations_names_dict,
                                    replicate_1_column_names,
                                    replicate_2_column_names,
                                    title_fontsize,
                                    x_dimension = 10,
                                    do_print = False,
                                    do_scale = False
                                    ):
    for index in index_selection:
        rna_counts_rep_1, rna_counts_rep_2, present_mutations_order_rna = \
            perturbations.subset_counts_five_perturbations_v2(
                        mutation_summary_df.loc[index],
                        rna_counts_df,
                        perturbations_order,
                        replicate_1_column_names,
                        replicate_2_column_names)
        gdna_counts_rep_1, gdna_counts_rep_2, present_mutations_order_gDNA = \
            perturbations.subset_counts_five_perturbations_v2(
                        mutation_summary_df.loc[index],
                        gdna_counts_df,
                        perturbations_order,
                        replicate_1_column_names,
                        replicate_2_column_names)

        gDNA_tick_labels = [perturbations_names_dict[x] for x in present_mutations_order_gDNA]
        visualize_four_mpra_arrays(
                                 gdna_counts_rep_1.to_numpy(),
                                 gdna_counts_rep_2.to_numpy(),
                                 rna_counts_rep_1.to_numpy(),
                                 rna_counts_rep_2.to_numpy(),
                                 gdna_counts_rep_1.shape[0],
                                 ytick_labels = gDNA_tick_labels,
                                 x_dimension = x_dimension,
                                 do_scale=do_scale,
                                 do_return=True,
                                 label = "%s\n%s" %
                                         (mutation_summary_df.loc[index, 'gene_name'],
                                          mutation_summary_df.loc[index, 'fr_name']),
                                 title_fontsize = title_fontsize)
    if do_print:
        sub_corr_df = mutation_summary_df.loc[index_selection]
        for index, row in sub_corr_df.iterrows():
            current_string = "%s\t%s\t%s" % (row['gene_name'], row['fr_name'], row['source'])
            print(current_string)


def visualize_mutations_one_by_one(chosen_mutations_df,
                                   figsize=(9, 3),
                                   s=7,
                                   alpha=0.5,
                                   ):
    indexed_df = chosen_mutations_df.copy()
    indexed_df = indexed_df.reset_index()
    for i in indexed_df.index:
        cur_df = indexed_df.loc[i]
        visualize_mutation_prob(cur_df,
                               figsize = figsize,
                               name="%s %s" % (cur_df['gene_name'], cur_df['changed_loop']),
                               s = s, alpha = alpha)


def print_chosen_mutations(mutation_choice_df):
    for index, row in mutation_choice_df.iterrows():
        print("Gene: ", row['gene_name'])
        print("Source: ", row['source'])
        print("Fragment name: ", row['fr_name'])
        print("Loop: ", row['changed_loop'])
        original_sequence = row['original']
        mutated_sequence = row['mutated']
        rescued_sequence = row['rescued']
        first_loop_string = row['first_loop'].get_dot_bracket_string(np.zeros(len(original_sequence)))
        second_loop_string = row['second_loop'].get_dot_bracket_string(np.zeros(len(original_sequence)))
        mutated_sequence_difference = "".join(
            ["." if x == y else x for x, y in zip(mutated_sequence, original_sequence)])
        rescued_sequence_difference = "".join(
            ["." if x == y else x for x, y in zip(rescued_sequence, original_sequence)])

        print("Sequence")
        print(original_sequence)
        print("First loop")
        print(first_loop_string)
        print("Second loop")
        print(second_loop_string)
        print("Mutated sequence difference")
        print(mutated_sequence_difference)
        print(mutated_sequence)
        print("Rescued sequence difference")
        print(rescued_sequence_difference)
        print(rescued_sequence)
        print('\n')


def visualize_pert_profile_single_array(
                                    index_selection,
                                    corr_df,
                                    counts_df,
                                    perturbations_order,
                                    column_names,
                                    cmap = 'YlOrRd',
                                    do_print = False,
                                    do_scale = False,
                                    do_return = False
                                    ):
    for index in index_selection:
        counts = perturbations.subset_counts_five_perturbations_from_averaged(
                        corr_df.loc[index],
                        counts_df,
                        perturbations_order,
                        column_names)
        counts = counts.to_numpy()
        if do_scale:
            for i in range(counts.shape[0]):
                #counts[i] = (counts[i] - counts[i].min()) / (counts[i].max() - counts[i].min())
                counts[i] = counts[i] / counts[i].sum()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 3))
        ax.imshow(counts,
                  cmap=cmap)
        #ax.set_yticks(np.arange(8))
        if do_return:
            return ax
        plt.show()
        # if not yticks_labels is None:
        #     ax.set_yticklabels(yticks_labels, fontdict=None, minor=False)
        # ax.set_title("Replicate 1: \n %s" % label, fontsize=title_fontsize)


    if do_print:
        sub_corr_df = corr_df.loc[index_selection]
        for index, row in sub_corr_df.iterrows():
            current_string = "%s\t%s\t%s" % (row['gene_name'], row['fr_name'], row['source'])
            print(current_string)


def visualize_pert_profile_single_array_separately(
                                    index_selection,
                                    corr_df,
                                    counts_df,
                                    perturbations_order,
                                    column_names,
                                    cmap = 'YlOrRd',
                                    do_print = False,
                                    do_scale = False,
                                    do_return = False
                                    ):
    for index in index_selection:
        counts = perturbations.subset_counts_five_perturbations_from_averaged(
                        corr_df.loc[index],
                        counts_df,
                        perturbations_order,
                        column_names)
        counts = counts.to_numpy()
        if do_scale:
            for i in range(counts.shape[0]):
                #counts[i] = (counts[i] - counts[i].min()) / (counts[i].max() - counts[i].min())
                counts[i] = counts[i] / counts[i].sum()

        fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(5, 3))
        vmin = np.quantile(counts, 0.001)
        vmax = np.quantile(counts, 0.999)
        print("Vmin: %.3f, vmax: %.3f" % (vmin, vmax))
        for i in range(len(perturbations_order)):
            axs[i].imshow(counts[i].reshape(counts.shape[1], 1).transpose(),
                  cmap=cmap, vmin = vmin, vmax = vmax)
        #ax.set_yticks(np.arange(8))
        if do_return:
            return axs
        plt.show()
        # if not yticks_labels is None:
        #     ax.set_yticklabels(yticks_labels, fontdict=None, minor=False)
        # ax.set_title("Replicate 1: \n %s" % label, fontsize=title_fontsize)


    if do_print:
        sub_corr_df = corr_df.loc[index_selection]
        for index, row in sub_corr_df.iterrows():
            current_string = "%s\t%s\t%s" % (row['gene_name'], row['fr_name'], row['source'])
            print(current_string)

