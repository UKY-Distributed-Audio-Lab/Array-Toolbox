
[r, c, n] = size(gout);
ulimim = .4;
for k=1:n
      gdum = squeeze(gout(:,:,k));
       thr = abs(min(min(gdum)));
       if thr > ulimim;
           gdum = zeros(size(gdum));
           thr = 0;
       else
        gdum = gdum-3*thr;
       end
       
      imagesc(fovxa, fovya, gdum, [0 ulimim])
      colormap(gray)
      axis('xy')
      xlabel('meters')
      ylabel('meters')
%       set(gcf,'Position', [232+10   258+10   764   420])
%       hold on
%       plot(envpar.mpos(1,:),envpar.mpos(2,:),'kv', 'MarkerSize', 14) 
%       plot(fovxa(xm),fovya(ym),'co','MarkerSize',14)
%       plot(imposesx(kv),imposesy(kv),'cx','MarkerSize',14)
%       hold off
      mmv(k) = getframe;
end