/** This file contains settings for the simulation that users can modify if they have the need.
 * To apply the modified settings the simulation tool must be re-compiled.
 * To modify a setting, replace the number (eg 1000) with the new value without modifying the name
 */
 
/** The largest contiguous block of memory that could be requested in the
 * process of simulation is CONTIG_WIDTH integers. This setting's default value 
 * is 1000.
 *
 * This could be decreased to help long simulations or lower-end machines. 
 * Increasing this may provide some speed gain. 
 */
#define CONTIG_WIDTH 1000
 
 /** The maximum number of characters allowed in a name field.
 * These include names of SNPs, names of genotypes loaded from files, 
 * names of generated genotypes, and save-as-you-go filenames. Default is 30.
 *
 * Increase this if there is a risk some names may be longer than 
 * this value.
 */
#define NAME_LENGTH 30
 
