function date = sessionToDate(session)
    date = split(session, 'd');
    date = date{2};
    year = str2double(date(1:4));
    month = str2double(date(5:6));
    day = str2double(date(7:8));
    date = datetime(year, month, day);
end