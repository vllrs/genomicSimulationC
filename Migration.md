Migration Guide        {#migration}
============================

Sometimes, function signatures may need to be changed at a version update. The changes needed for a script to be updated to run in a newer version of genomicSimulation can be found in this guide. Angled brackets represent function inputs.

This guide belongs to the C version of genomicSimulation. Function signature changes in the R version of genomicSimulation will be managed with `.Deprecated` warnings.

# 0.2.5 → 0.2.6

Function signatures for file output functions were simplified. So were function signatures of breeding value and allele count matrix calculation functions. One extra parameter was added to `load_data_files`.

<table>
<tr><th>0.2.5<th>0.2.6
<tr>
	<td>`load_data_files(<d>, <data_file>, <map_file>, <effect_file>)`
	<td>`load_data_files(<d>, <data_file>, <map_file>, <effect_file>, DETECT_FILE_FORMAT)`
<tr>
	<td>`save_markerblocks(<f>, <d>, <b>)` 
	<td>`save_markerblocks(<!fname>, <d>, <b>, NO_MAP)`<br>followed by running the command ```sed 's/^/0\t0\tb0\tb\t/' <!fname> | sed '1i\Chrom\tPos\tName\tClass\tMarkers'```<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)<br>(It is suggested you use one of the simpler new output formats without meaningless columns instead.)
<tr>
	<td>`save_group_genotypes(<f>, <d>, <group>)`
	<td>`save_genotypes(<!fname>, <d>, <group>, GSC_FALSE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_transposed_group_genotypes(<f>, <d>, <group>)`
	<td>`save_genotypes(<!fname>, <d>, <group>, GSC_TRUE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_count_matrix(<f>, <d>, <allele>)`
	<td>`save_allele_counts(<!fname>, <d>, NO_GROUP, <allele>, GSC_FALSE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_group_count_matrix(<f>, <d>, <allele>, <group>)`
	<td>`save_allele_counts(<!fname>, <d>, <group>, <allele>, GSC_FALSE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_one_step_pedigree(<f>, <d>)`
	<td>`save_pedigrees(<!fname>, <d>, NO_GROUP, GSC_FALSE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_group_one_step_pedigree(<f>, <d>, <group>)`
	<td>`save_pedigrees(<!fname>, <d>, <group>, GSC_FALSE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_full_pedigree(<f>, <d>, <group>)`
	<td>`save_pedigrees(<!fname>, <d>, NO_GROUP, GSC_TRUE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_group_full_pedigree(<f>, <d>, <group>)`
	<td>`save_pedigrees(<!fname>, <d>, <group>, GSC_TRUE)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)

<tr>
	<td>`save_bvs(<f>, <d>, <effID>)`
	<td>`save_bvs(<!fname>, <d>, NO_GROUP, <effID>)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`save_group_bvs(<f>, <d>, <group>, <effID>)`
	<td>`save_bvs(<!fname>, <d>, <group>, <effID>)`<br>(where `<!fname>` is the file name that would, in the old version, have been opened as file pointer `<f>`)
<tr>
	<td>`calculate_bvs(<d>->m, <eff>)`<br>where `<d>` is a SimData object
	<td>`calculate_bvs(<d>, NO_GROUP, <eff>)`
<tr>
	<td>`calculate_bvs(<m>, <eff>)`<br>
	if `<m>` does not belong to a SimData object
	<td>`BidirectionalIterator it = gsc_create_bidirectional_iter_fromAM(<m>, NO_GROUP);`<br>`gsc_calculate_utility_bvs(&it, <eff>);`<br>`gsc_delete_bidirectional_iter(&it);`
<tr>
	<td>`calculate_group_bvs(<d>, <group>, <eff>)`
	<td>`calculate_bvs(<d>, <group>, <eff>)`
<tr>
	<td>`calculate_count_matrix(<d>->m, <allele>, <counts>)`<br>where `<d>` is a SimData object
	<td>`gsc_delete_dmatrix(&<counts>);` (if needed)<br>`<counts> = calculate_allele_counts(<d>, NO_GROUP, <allele>);`
<tr>
	<td>`calculate_group_count_matrix(<d>->m, <group>, <allele>, &<counts>)`<br>where `<d>` is a SimData object
	<td>`gsc_delete_dmatrix(&<counts>);` (if needed)<br>`<counts> = calculate_allele_counts(<d>, <group>, <allele>);`
<tr>
	<td>`calculate_full_count_matrix(<d>->m, <allele>)`<br>where `<d>` is a SimData object
	<td>`calculate_allele_counts(<d>, NO_GROUP, <allele>)`
<tr>
	<td>`calculate_full_count_matrix(<m>, <allele>)`<br>if `<m>` does not belong to a SimData object
	<td>`DecimalMatrix counts = gsc_generate_zero_dmatrix(<m>->n_genotypes, <m>->n_markers);`<br>`gsc_calculate_utility_allele_counts(<m>->n_markers, <m>->n_genotypes, <m>->alleles, <allele>, &counts);`<br>The result is in `counts`
<tr>
	<td>`gsc_save_names_header(<f>, <n>, <names>)`
	<td>`gsc_save_utility_genotypes(<f>, NULL, <n>, <names>, GSC_FALSE)`
<tr>
	<td>`gsc_save_allele_matrix(<f>, <m>)`
	<td>`BidirectionalIterator it = gsc_create_bidirectional_iter_fromAM(<m>, NO_GROUP);`<br>`gsc_save_utility_genotypes(<f>, &it, 0, NULL, GSC_FALSE);`<br>`gsc_delete_bidirectional_iter(&it);`
<tr>
	<td>`gsc_save_transposed_allele_matrix(<f>, <m>, <names>)`
	<td>`BidirectionalIterator it = gsc_create_bidirectional_iter_fromAM(<m>, NO_GROUP);`<br> `gsc_save_utility_genotypes(<f>, &it, <m>->n_markers, <names>, GSC_TRUE);`<br>`gsc_delete_bidirectional_iter(&it);`
<tr>
	<td>`gsc_save_allelematrix_full_pedigree(<f>, <m>, <parents>)`
	<td>`BidirectionalIterator it = gsc_create_bidirectional_iter_fromAM(<m>, NO_GROUP);`<br> `gsc_save_utility_pedigrees(<f>, &it, GSC_TRUE, <parents>->m);`<br>`gsc_delete_bidirectional_iter(&it);` 
<tr>
	<td>`gsc_DecimalMatrix bvs = gsc_calculate_bvs(<m>, <eff>);`<br>`gsc_save_manual_bvs(<f>, &bvs, <m>->ids, (const char**) <m>->names);`<br>`gsc_delete_dmatrix(&bvs);`
	<td>`gsc_BidirectionalIterator it = gsc_create_bidirectional_iter_fromAM(<m>, NO_GROUP);`<br>`gsc_save_utility_bvs(<f>, &it, <eff>);`<br>`gsc_delete_bidirectional_iter(&it);`
</table>


# 0.2.4 → 0.2.5

See the file `names_stopgap.h` in the v0.2.5 release.

```
#define GSC_DEPRECATED_VERBOSE
#include "names_stopgap.h"
```
