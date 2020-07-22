function Ya = create_Ya2(Y,id,a)

Ya{a}{1} = squeeze(Y(:, a, :,1)); 
Ya{a}{1} = Ya{a}{1}(id,:);
           
Ya{a}{2} = squeeze(Y(:, a, :,2)); 
Ya{a}{2} = Ya{a}{2}(id,:);
            
Ya{a}{3} = squeeze(Y(:, a, :,3)); 
Ya{a}{3} = Ya{a}{3}(id,:);
            
Ya{a}{4} = squeeze(Y(:, a, :,4)); 
Ya{a}{4} = Ya{a}{4}(id,:);

end
