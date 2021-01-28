function result = getPNStdevBhandawat(PNs)
            % for a numPNs x numOdors matrix, returns the estimated
            % standard deviation for each response, based on Bhandawat et
            % al 2007 Fig 1e.
            
            numPNs = size(PNs,1);
            numOdors = size(PNs,2);
            result = zeros(numPNs, numOdors);
            
            % column 1 is bin centers (mean spike rate in Hz)
            % column 2 is the average st.dev. for 50 ms windows with that
            % mean spike rate (in Hz)
            BhandawatData = [
                10	2.929936306
                30	6.904458599
                50	8.687898089
                70	10.31847134
                90	11.2611465
                110	11.69426752
                130	11.69426752
                150	10.70063694
                170	9.78343949
                190	9.732484076
                210	8.866242038
                230	8.687898089
                250	7.363057325
                270	8.152866242
                290	10.67515924
                310	9.910828025];
            
            for i=1:numPNs
                for j=1:numOdors
                    % find index of the bin center closest to this PN
                    % response value
                    [~, index] = min(abs(BhandawatData(:,1)-PNs(i,j)));
                    result(i,j) = BhandawatData(index,2);
                end
            end
            
            
end
