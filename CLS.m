load('Z:\Users\John\rhenium complex\2Dspecs\ReCN_THF_absorptive_hf\AbsorptiveSpecs.mat') 

settings.doBatch = 0;
settings.doSingle = 1;
settings.plotAbsorptive = 0;
settings.findRange = 1;

N = 13;
AbsFinal = output.Absorptive(:,:,N);
AbsFinal2=output.Absorptive(:,:,N);
wn = output.wn;
w3wn = output.w3wn;

AbsFinal2 = flipud(AbsFinal2);
AbsFinal2 = fliplr(AbsFinal2);

xrange = 106:111;
yrange = 165:145

if settings.plotAbsorptive
    for n = 1:size(AbsFinal,3)
        numSteps = 25;
        figure(12)
        contourf(wn,w3wn,AbsFinal(:,:,n),numSteps)

        xlim([2000 2040]);
        ylim([2000 2040]);
        axis square 
    end
end

if settings.findRange
    figure(13)
    contourf(AbsFinal2(:,:), 20)
end


[C,I] = min(AbsFinal);
I = squeeze(I);



if settings.doSingle
    [Islope, inter] = polyfit(xrange,I(xrange),1);
    slope = Islope(1,1);
    display(slope) 
    figure(10) 
    plot(xrange,I(xrange))
end

if settings.doBatch         %does not work
    for n=1:size(I,2)
        [Islope, inter] = polyfit(xrange, I(range), n);
        slope = Islope(:,n);
        slopeAll(n)=slope;
        slopeAll = slopeAll';
    end
end

