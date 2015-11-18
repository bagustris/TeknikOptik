% image correction
raw=aviread(dr(1).name,1);   % read first image
bg=double(raw.cdata(:,:,1)); % rgb all same for grayscale movie
simage([im bg]);  % display image in false color
title('Image of an Animal to be Tracked and Background.');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;
