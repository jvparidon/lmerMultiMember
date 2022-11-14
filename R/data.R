#' National Football League scores from the 2021 season
#'
#' A subset of NFL scores from the FiveThirtyEight's football model dataset,
#' released under CC-BY license and reproduced here with attribution.
#'
#' @format ## `nfl_scores_2021`
#' A data frame with 285 rows and 7 columns:
#' \describe{
#'   \item{date}{Date on which a game was played}
#'   \item{season}{NFL season in which a game was played}
#'   \item{home_team, visiting_team}{3-letter codes for the playing teams}
#'   \item{home_score, visiting_score}{Scores for the playing teams}
#'   \item{winner}{Whether the home team or the visiting team won the game}
#' }
#' @source <https://github.com/fivethirtyeight/data/tree/master/nfl-elo>
"nfl_scores_2021"


#' ATP men's doubles tennis matches from the 2010 to 2019
#'
#' A subset of ATP men's doubles tennis data compiled by Jeff Sackmann,
#' released under CC-BY-NC-SA license and reproduced here with attribution.
#'
#' @format ## `atp_doubles`
#' A data frame with 13038 rows and 47 columns:
#' \describe{
#'   \item{tourney_id}{Unique identifier for the tournamens}
#'   \item{tourney_name}{Name of the tournament}
#'   \item{surface}{Surface the match was played on}
#'   \item{draw_size}{???}
#'   \item{tourney_level}{???}
#'   \item{tourney_date}{Date on which the tournament started}
#'   \item{match_num}{???}
#'   \item{score}{Scores per set}
#'   \item{best_of}{Maximum number of sets, either best-of-3 or best-of-5}
#'   \item{round}{Round in the tournament in which the match was played}
#'   \item{team1_win}{1 if team 1 won, 0 if team 2 won}
#'   \item{team1_player1_id,
#'         team1_player2_id,
#'         team2_player1_id,
#'         team2_player2_id}{Unique identifier for each player on each team}
#'   \item{team1_seed, team2_seed}{Position at which each team was seeded}
#'   \item{team1_entry, team2_entry}{???}
#'   \item{team1_player1_name,
#'         team1_player2_name,
#'         team2_player1_name,
#'         team2_player2_name}{Name for each player on each team}
#'   \item{team1_player1_hand,
#'         team1_player2_hand,
#'         team2_player1_hand,
#'         team2_player2_hand}{Handedness for each player on each team}
#'   \item{team1_player1_ht,
#'         team1_player2_ht,
#'         team2_player1_ht,
#'         team2_player2_ht}{Height in cm for each player on each team}
#'   \item{team1_player1_ioc,
#'         team1_player2_ioc,
#'         team2_player1_ioc,
#'         team2_player2_ioc}{IOC for which each player on each team competes}
#'   \item{team1_player1_age,
#'         team1_player2_age,
#'         team2_player1_age,
#'         team2_player2_age}{Age for each player on each team}
#'   \item{team1_player1_rank,
#'         team1_player2_rank,
#'         team2_player1_rank,
#'         team2_player2_rank}{ATP ranking for each player on each team}
#'   \item{team1_player1_rank_points,
#'         team1_player2_rank_points,
#'         team2_player1_rank_points,
#'         team2_player2_rank_points}{???}
#' }
#' @source <https://github.com/JeffSackmann/tennis_atp>
"atp_doubles"
