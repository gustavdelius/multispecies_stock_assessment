local has_mermaid = false

local function escape_html(s)
  return (s:gsub("&", "&amp;")
           :gsub("<", "&lt;")
           :gsub(">", "&gt;"))
end

function CodeBlock(el)
  if el.classes:includes("mermaid") and quarto.doc.is_format("html:js") then
    has_mermaid = true
    local html = string.format(
      '<div class="cell" data-layout-align="default">\n' ..
      '<div class="cell-output-display">\n' ..
      '<div>\n' ..
      '<figure class="figure">\n' ..
      '<div>\n' ..
      '<pre class="mermaid mermaid-js">%s</pre>\n' ..
      '</div>\n' ..
      '</figure>\n' ..
      '</div>\n' ..
      '</div>\n' ..
      '</div>',
      escape_html(el.text)
    )
    return pandoc.RawBlock("html", html)
  end
end

function Meta(meta)
  if has_mermaid and quarto.doc.is_format("html:js") then
    local share_path = os.getenv("QUARTO_SHARE_PATH")
    if share_path then
      local mermaid_dir = pandoc.path.join({share_path, "formats", "html", "mermaid"})
      quarto.doc.add_html_dependency({
        name = "quarto-diagram",
        scripts = {
          {
            name = "mermaid.min.js",
            path = pandoc.path.join({mermaid_dir, "mermaid.min.js"})
          },
          {
            name = "mermaid-init.js",
            path = pandoc.path.join({mermaid_dir, "mermaid-init.js"}),
            afterBody = true
          }
        },
        stylesheets = {
          {
            name = "mermaid.css",
            path = pandoc.path.join({mermaid_dir, "mermaid.css"})
          }
        }
      })
    end
  end
end

return {
  { CodeBlock = CodeBlock },
  { Meta = Meta }
}
