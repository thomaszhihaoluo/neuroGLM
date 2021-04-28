function Cells = prepare_Cells_for_neuroGLM(Cells, Sessions)

if iscell(Cells.rat) || isstring(Cells.rat)
    Cells.rat = Cells.rat{1};
end
Cells.sessid = Cells.sessid(1);
Cells.sess_date = datestr(Sessions.sess_date, 'yyyy_mm_dd');
Cells.sess_time = [];
Cells.nTrials = Sessions.nTrials;

Trials = PB_make_Trials(Sessions);
trial_index = num2cell(1:numel(Trials.leftBups))';
t_clicks_on = num2cell(Trials.stateTimes.clicks_on);

Trials.stateTimes.left_click_times = cell2mat(Trials.leftBups')';
n_left = cellfun(@numel, Trials.leftBups, 'uni', 0);
left_click_trial = cellfun(@(x,y) repmat(x,y,1), trial_index, n_left, 'uni', 0);
Trials.stateTimes.left_click_trial = cell2mat(left_click_trial);
Trials.stateTimes.left_clicks = cell2mat(cellfun(@(x,y) x+y, Trials.leftBups, t_clicks_on, 'uni', 0)')';

Trials.stateTimes.right_click_times = cell2mat(Trials.rightBups')';
n_left = cellfun(@numel, Trials.rightBups, 'uni', 0);
right_click_trial = cellfun(@(x,y) repmat(x,y,1), trial_index, n_left, 'uni', 0);
Trials.stateTimes.right_click_trial = cell2mat(right_click_trial);
Trials.stateTimes.right_clicks = cell2mat(cellfun(@(x,y) x+y, Trials.rightBups, t_clicks_on, 'uni', 0)')';

Cells.Trials = Trials;