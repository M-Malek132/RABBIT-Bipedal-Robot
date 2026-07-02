function update_rabbit_frame(q,handles)

[stance_foot, swing_foot, hip, ...
 stance_knee, swing_knee, torso_top] = rabbit_kinematics(q);

% stance leg
set(handles.stance_leg,'XData',[hip(1) stance_knee(1) stance_foot(1)], ...
                       'YData',[hip(2) stance_knee(2) stance_foot(2)]);

% swing leg
set(handles.swing_leg,'XData',[hip(1) swing_knee(1) swing_foot(1)], ...
                      'YData',[hip(2) swing_knee(2) swing_foot(2)]);

% torso
set(handles.torso,'XData',[hip(1) torso_top(1)], ...
                  'YData',[hip(2) torso_top(2)]);

% joints
set(handles.hip,'XData',hip(1),'YData',hip(2));
set(handles.stance_knee,'XData',stance_knee(1),'YData',stance_knee(2));
set(handles.swing_knee,'XData',swing_knee(1),'YData',swing_knee(2));
set(handles.stance_foot,'XData',stance_foot(1),'YData',stance_foot(2));
set(handles.swing_foot,'XData',swing_foot(1),'YData',swing_foot(2));

end
