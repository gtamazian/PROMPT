function hinges_found = pdb_extract_hinges(pdb_begin_filename, ...
    pdb_end_filename, begin_conf_numb, end_conf_numb, out_path, ...
    hinge_angle_threshold)

[hinges_ang_ranges, hinges_types] = find_hinges(pdb_begin_filename, ...
    pdb_end_filename, begin_conf_numb, end_conf_numb, hinge_angle_threshold);

if isempty(hinges_ang_ranges)
    hinges_found = false;
else
    extract_hinges(pdb_begin_filename, pdb_end_filename, begin_conf_numb, ...
        end_conf_numb, out_path, hinges_ang_ranges, hinges_types);
    hinges_found = true;
end
end

function [hinges, htypes] = find_hinges(pdb_beg_filename, ...
    pdb_end_filename, beg_conf_numb, end_conf_numb, hinge_angle_threshold)

ANGLE_THRESHOLD = 0.5;
NEIGHBOUR_THRESHOLD = 30;
TAIL_THRESHOLD = 30;

% FINDHINGES finds protein regions which are most changeable during protein
% transformation between two conformations.
% It returns list of pairs of all hinges in protein. Element of
% a pair is a begining and ending of hinge. They are ordering numbers of
% atoms of main chain from protein n-terminus.

if ~exist('hinge_angle_threshold', 'var') || ~hinge_angle_threshold
    hinge_angle_threshold = ANGLE_THRESHOLD;
end

pdb_struct_beg = pdbread(pdb_beg_filename);
pdb_struct_end = pdbread(pdb_end_filename);

first_conf = pdb_struct_beg;
first_conf.Model = first_conf.Model(beg_conf_numb);
last_conf = pdb_struct_end;
last_conf.Model = last_conf.Model(end_conf_numb);

trm = trmcreate(first_conf, last_conf, 0);

angdiff = abs(circdist(trm.psi(:,1), trm.psi(:,end)));
ind = 1:length(angdiff);

hinges = [];
htypes = [];

candidates = ind(angdiff > hinge_angle_threshold);

if ~isempty(candidates)
    neighbours = candidates(2:end) - candidates(1:end-1);
    ind = 1:length(neighbours);
    borders = ind(neighbours > NEIGHBOUR_THRESHOLD);
    lowind = [1 borders + 1];
    highind = [borders length(candidates)];
    hinges = [candidates(lowind).' candidates(highind).'];
    hinges = hinges(hinges(:,1) ~= hinges(:,2),:);

    htypes = repmat({'loop'}, size(hinges,1), 1);
    if hinges(1,1) <= TAIL_THRESHOLD
        htypes{1} = 'n-tail';
        hinges(1,1) = 1;
    end
    if hinges(end,end) >= size(trm.psi,1) - TAIL_THRESHOLD
        htypes{end} = 'c-tail';
        hinges(end,end) = size(trm.psi,1);
    end
end
end

function hinge_filenames = extract_hinges(first_conf_pdb_path, ...
    last_conf_pdb_path, first_conf_ind, last_conf_ind, out_path, ...
    hinges, hinge_type)

first_conf = pdbread(first_conf_pdb_path);
first_conf.Model = first_conf.Model(first_conf_ind);
first_conf = pdbbackbone(first_conf);

last_conf = pdbread(last_conf_pdb_path);
last_conf.Model = last_conf.Model(last_conf_ind);
last_conf = pdbbackbone(last_conf);

hinges_count = size(hinges, 1);
indices = (1:hinges_count).';
hinge_filenames = repmat({''}, hinges_count, 1);
begin_conf_pdb_path =repmat(first_conf_pdb_path, hinges_count, 1);
end_conf_pdb_path = repmat(last_conf_pdb_path, hinges_count, 1);
begin_conf_num = repmat(first_conf_ind, hinges_count, 1);
end_conf_num = repmat(last_conf_ind, hinges_count, 1);

for i = 1:size(hinges, 1)
    hinge_first_conf = ...
        pdbgetatomsrange(first_conf,hinges(i,1),hinges(i,end)+3);
    hinge_last_conf = ...
        pdbgetatomsrange(last_conf,hinges(i,1),hinges(i,end)+3);
    trm = trmcreate(hinge_first_conf, hinge_last_conf, 0, [], []);
    trm.U = trmsuperpos(trm);
    pdb = trm2pdb(trm, hinge_first_conf);
    filename = ['hinge_' num2str(i)];
    hinge_filenames{i} = fullfile(out_path,[filename '.pdb']);
    pdbwrite(hinge_filenames{i}, pdb);
end

hinge_n_tail_torsion_angle = hinges(:,1);
hinge_c_tail_torsion_angle = hinges(:,2);
T = table(indices,hinge_n_tail_torsion_angle, ...
    hinge_c_tail_torsion_angle, ...
    hinge_filenames, begin_conf_pdb_path, ...
    end_conf_pdb_path, begin_conf_num, ...
    end_conf_num, hinge_type);
writetable(T, fullfile(out_path,'hinges.csv'));

end

function pdb = pdbgetatomsrange(pdb, begin_atom, end_atom)
    pdb.Model.Atom = pdb.Model.Atom(begin_atom:end_atom);
end
