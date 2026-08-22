function ch3_logln(logfile, msg)
%CH3_LOGLN  Print a march progress line, and append it to a log immediately.
%
%   ch3_logln(logfile, msg)
%   ch3_logln('', msg)          print only
%
% WHY A MARCH NEEDS A FILE AND NOT JUST fprintf.  MATLAB buffers stdout under
% -batch, so a campaign that runs for an hour prints nothing until it exits and
% is indistinguishable from a hang -- there is no way to tell a stalled solve
% from a slow one, or to see which rung it reached before it died.  The log is
% opened, written and CLOSED on every call, which is what forces the write
% through and makes the file readable live from another shell:
%
%       tail -f Results/ch3_leantall.log
%
% The timestamp goes to the FILE only.  Stdout gets the bare message, because
% when you are watching interactively the timestamps are noise, and when you
% are reading the file afterwards they are the whole point -- they turn the log
% into per-stage timings.
%
% A FAILED LOG MUST NOT KILL THE MARCH.  This is the same argument ch3_col_solve
% makes for its checkpoint OutputFcn, and it applies with more force here: a
% march is hours of solves, the log path routinely comes from a p loaded out of
% a .mat written on another machine, and erroring at rung six over a stale
% scratch directory throws away everything the run had done.  So a failure to
% open warns ONCE and the march continues -- the message still reaches stdout,
% so nothing is lost but the file copy.
%
% NOTE this differs from the two predecessor helpers it replaces, which both
% raised an error.  In practice they failed at their first call, before any
% solve, so the fail-fast cost nothing; the risk is only realized when the path
% goes bad mid-run, which is exactly when erroring is most expensive.
%
% Inputs
%   logfile : path to append to; '' or [] prints without writing
%   msg     : the line (or multi-line block) to record.  A trailing newline is
%             optional -- callers that build blocks with sprintf may or may not
%             leave one, and both read the same in the log.
%
% See also CH3_LEAN_TALL_MARCH, CH3_REALIZABILITY_MARCH, CH3_IMPACT_MARCH.

msg = char(msg);
msg = regexprep(msg, '\n+$', '');       % callers differ on the trailing newline

fprintf('%s\n', msg);

if nargin < 1 || isempty(logfile)
    return;
end

persistent warned
fid = fopen(logfile, 'a');
if fid < 0
    if isempty(warned)
        warning('ch3_logln:open', ...
                ['Cannot open the log file "%s" for append; continuing with ' ...
                 'console output only.\nSet the log path to somewhere ' ...
                 'writable, or '''' to disable the file copy.'], logfile);
        warned = true;
    end
    return;
end

fprintf(fid, '%s  %s\n', datestr(now, 'HH:MM:SS'), msg); %#ok<TNOW1,DATST>
fclose(fid);            % close every time -- that is what makes it readable live

end
