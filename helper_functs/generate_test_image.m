function [im,ind,subInd] = generate_test_image(form,N1,N2)
% generates image im of size N1 x N2 with the requested form, as well as
% linear indices ind and subscript indices subInd on which to
% perform analysis
%
% currently only supports form = 'sqCirc'

ind = zeros(4,1);
savefolder = "Emp_Bayes_1D_2D/helper_functs/test_images/";
mkdir(savefolder);

% N_s is the length of the longer side, N_b the length of the shorter side
if N1 <= N2
    N_s = N1;
    N_b = N2;
else
    N_s = N2;
    N_b = N1;
end

% 'sqCirc' creates an image of a centered circle with a square inset
if strcmp(form,'sqCirc')
    h_base = 1;
    h_circ = 5;
    h_square = 3;
    r_circ = round(N_s/2-N_s/9);
    r_square = round(N_s/2-N_s/3);

    im = h_base*ones(N1,N2);
    center = [round((N1+1)/2) round((N2+1)/2)];
    for ii = 1:N1
        for jj = 1:N2
            if sqrt((ii-center(1))^2+(jj-center(2))^2) <= r_circ
                im(ii,jj) = h_circ;
            end
            if abs(ii-center(1)) <= r_square && abs(jj-center(2)) <= r_square
                im(ii,jj) = h_square;
            end
        end
    end
    im = im./max(im,[],'all');
    subInd = [round(N_s/2) round(N_s/9);
              round(N_s/3) round(N_s/2);
              round(N_s/2) round(N_s/2);
              3 3];
    ind = sub2ind([N1 N2],subInd(:,1),subInd(:,2));

    figure;imagesc(im);colorbar;
    savefile = sprintf("2d_testImage_sqCirc_%d_%d",N1,N2);
    save(strcat(savefolder,savefile),'im');
end