% * Analytic steps (simplified, least amount needed for paper):
%     1. Obtain functional activity for each electrode (nodal measures)
%         * TF power
%         *  set basic statistical threshold p < 0.05)
%         * Extract time-frequency clusters (contiguous in TF space)
%     2. Obtain anatomical information for each electrode
%         * standardized atlas-coordinate
%         * atlas-based ROI
%         * proximity to white-matter
%     3. Describe networks:
%         * Do additional functional testing for each activity pattern: e.g., left vs. right selectivity; instruct vs. wait selectivity; correlation with reaction time
%         * Then each activity pattern will be tagged with functional selectivity
%         * Label networks in terms of functional significance  (e.g., motor effector, response inhibition, salience)
%         * Describe anatomical pattern of each network
%     4. Study connectivity between networks
%         * Granger causality vs. phase synchrony to describe connectivity between functional networks





% * Code plan:
%     * Create new analysis folder motor_networks:  le_mn_..
%     * le_mn_spectralChange
%         * Goal: Obtain functional activity for each electrode (nodal measures; step 1)
%         * Returns a structure with clusters of spectral change from each electrode
%         * Also tag each cluster with subj, sess, and electrode (possibly,  anatomical information)
%         * Generalize to allow for tilt regression and phase reset
%         * Generalize to combine across electrodes for (ROI analysis or network analysis)
%     * le_mn_selectivity
%         * Takes a structure of cluster of activations from (spectralChange) and performs a variety of functional activations
%         * Returns the structure with results of functional battery additional functional
%         * So each activation pattern is tagged with functional information
%     * le_mn_clusterPatterns 
%         * Can cluster functional patterns based on spectral properties and functional selectivity
%     * le_mn_connectivity
%     * le_mn_getAnat
%         * Goal: gets anatomical data provided by Joel stein
%         * Can assign anatomical information to structure above
%     * le_mn_plotAnat
%         * Plotting networks on brain surface