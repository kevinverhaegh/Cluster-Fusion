pathdef='/media/DATA/GPT Simulations (all)/100partADK/';
names = cell(20,1);
names{1,1}='UCP001';
names{2,1}='UCP002';
names{3,1}='UCP003';
names{4,1}='UCP004';
names{5,1}='UCP005';
names{6,1}='UCP006';
names{7,1}='UCP007';
names{8,1}='UCP008';
names{9,1}='UCP009';
names{10,1}='UCP010';
names{11,1}='UCP011';
names{12,1}='UCP012';
names{13,1}='UCP013';
names{14,1}='UCP014';
names{15,1}='UCP015';
names{16,1}='UCP016';
names{17,1}='UCP017';
names{18,1}='UCP018';
names{19,1}='UCP019';
names{20,1}='UCP020';

for i=1:20
    names{i,1}=[pathdef, names{i,1}];
end

AutomizeSimsSort(names);